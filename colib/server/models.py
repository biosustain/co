from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy import func, CheckConstraint, select
from sqlalchemy.dialects import postgres
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.mutable import MutableDict
from sqlalchemy.orm import backref, deferred, aliased

db = SQLAlchemy()


component_removed_features_map = db.Table('component_removed_features_map',
    db.Column('component_id', db.Integer, db.ForeignKey('component.id'), primary_key=True),
    db.Column('feature_id', db.Integer, primary_key=True),
    db.Column('feature_component_id', db.Integer, primary_key=True),
    db.ForeignKeyConstraint(['feature_id', 'feature_component_id'], ['feature.id', 'feature.component_id']))


class Component(db.Model):
    """

    Feature Differences
    -------------------

    - Added features are simply added to the :class:`Component`; these are automatically combined with the inherited
      features from :attr:`parent`.
    - Changed features, excluding features where only `position` has changed, are represented by a deletion
      followed by an addition.
    - New positions for changed features are computed from the translation table `tt`. When the size of a feature has
      changed, it is expected that this change is represented by replacing the feature for the mutated component.

    """
    id = db.Column(db.Integer, primary_key=True)
    parent_id = db.Column(db.Integer, db.ForeignKey(id))

    sequence = deferred(db.Column(db.Text, nullable=True))

    #tt = db.Column(postgres.JSON)

    date_created = db.Column(db.DateTime(timezone=True), default=func.now(), nullable=False)

    # genbank_meta = db.Column(MutableDict.as_mutable(postgres.HSTORE))

    parent = db.relationship('Component', remote_side=[id], backref=backref('descendants', lazy='dynamic'))
    removed_features = db.relationship('BoundFeature', secondary=component_removed_features_map)

    def __init__(self, sequence=None, parent=None):
        self.sequence = sequence
        self.parent = parent

    @hybrid_property
    def is_mutated(self):
        return self.parent is not None

    @hybrid_property
    def size(self):
        return func.length(self.sequence)
    def get_features(self, include_inherited=True):
        if include_inherited:
            ancestors = db.session.query(Component.id, Component.parent_id)\
                .filter(Component.id==self.id)\
                .cte(name="ancestors", recursive=True)

            ancestors_alias = aliased(ancestors, name='a')
            components_alias = aliased(Component, name="c")
            ancestors = ancestors.union_all(
                db.session.query(components_alias.id, components_alias.parent_id)\
                    .filter(components_alias.id==ancestors_alias.c.parent_id))

            query = db.session.query(ancestors.c.id, ancestors.c.parent_id)\
             #   .group_by(ancestors.c.parent_id)

            return query
        else:
            return self.added_features
        # TODO map.

    def get_lineage(self, with_self=True):
        # WITH RECURSIVE included_components(parent_id, id) AS (
        #     SELECT parent_id, id FROM component WHERE id=2
        #   UNION ALL
        #     SELECT p.parent_id, p.id
        #     FROM included_components pr, component p
        #     WHERE p.id = pr.parent_id
        #   )
        #
        # SELECT parent_id, id
        # FROM included_components;
        ancestors = db.session.query(Component.id, Component.parent_id)\
            .filter(Component.id==self.id)\
            .cte(name="ancestors", recursive=True)

        ancestors_alias = aliased(ancestors, name='a')
        components_alias = aliased(Component, name="c")
        ancestors = ancestors.union_all(
            db.session.query(components_alias.id, components_alias.parent_id)
                .filter(components_alias.id==ancestors_alias.c.parent_id))

        filter_id = ancestors.c.id if with_self else ancestors.c.parent_id
        query = Component.query.filter(Component.id==filter_id)

        print()
        print(query)
        print()
        print()

        return query.all()

    def find_features(self, type=None, name=None, qualifiers=None, xref=None):
        pass

    def mutate(self, mutations):
        """


        """
        pass

    def __repr__(self):
        return "Component({})".format(self.id)


class BoundMutation(db.Model):
    """
    A :class:`BoundMutation` encodes all possible mutations on a component.

    .. note::

        The "Bound* in :class:`BoundMutation` and :class:`BoundFeature` is to indicate that these models are associated
        with a specific :class:`Component`. In contrast, a simple `Mutation` could refer to any sequence.

    .. attribute:: is_substitution

        `True` if both range and replacement are of length `1`.

    .. attribute:: is_insertion

        `True` if the replacement length is `> 1`.

        In an insertion, :attr:`start_offset` and :attr:`end_offset` are exclusive.

        For instance in insertion between 41..42, the replacement sequence is inserted after base 41.

    .. attribute:: is_deletion

        Always `True` if there is no replacement. If combined with :attr:`is_insertion`, it is interpreted
        as a deletion followed by an insertion.

    .. attribute:: is_duplication

        `True` if `new_sequence` is identical to two copies of the sequence present at the original range.

    .. attribute:: context

        Stores information about how this :class:`BoundMutation` was generated. For example:::

            `'{"range": {"start": ["FeatureStart", 1], "end": ["FeatureEnd", 1], "replacement": null}'`

        Context, if provided, should provide enough information to re-compute the :attr:`start`, :attr:`end` and
        :attr:`new_sequence` of the mutation.
    """
    id = db.Column(db.Integer, primary_key=True)
    component_id = db.Column(db.Integer, db.ForeignKey(Component.id), nullable=False)

    position = db.Column(db.Integer, nullable=False, index=True)
    size = db.Column(db.Integer, nullable=False, index=True)

    new_sequence = deferred(db.Column(db.Text, nullable=True))

    #context = db.Column(postgres.JSON)

    component = db.relationship(Component, backref=backref('mutations', lazy='dynamic'))

    __tablename__ = 'mutation'
    __table_args__ = (
        CheckConstraint(position >= 0),
        CheckConstraint(size >= 0),
    )

    @hybrid_property
    def orig_sequence(self):
        return func.substr(self.component.parent.sequence, self.start, self.end)

    @hybrid_property
    def new_size(self):
        return func.length(self.new_sequence)


class BoundFeature(db.Model):
    """
    :class:`BoundFeature` positions are *always* relative to :attr:`component`. Positions are zero-indexed.

    Linked components
    -----------------
    When another component is inserted into a :class:`Component`, a feature is created that links
    to the component. This link is stored in :attr:`component_link`.

    When a feature is changed, the link has to be flagged as broken. A broken link indicates that the feature sequence
    no longer matches the component sequence precisely.

    Caching
    -------
    Feature positions and sequences should be cached outside the database.
    """
    id = db.Column(db.Integer, primary_key=True)
    component_id = db.Column(db.Integer, db.ForeignKey(Component.id), primary_key=True)

    position = db.Column(db.Integer, nullable=False, index=True)
    size = db.Column(db.Integer, nullable=False, index=True)

    type = db.Column(db.String, index=True)  # Sequence Ontology terms preferred.
    name = db.Column(db.String, index=True)

    #meta = db.Column(postgres.JSON, nullable=False)

    source_link_id = db.Column(db.Integer, db.ForeignKey(Component.id))
    source_link_is_broken = db.Column(db.Boolean, default=False)

    component = db.relationship(Component,
                                foreign_keys=[component_id],
                                backref=backref('added_features', lazy='dynamic', order_by=position))

    source_link = db.relationship(Component, foreign_keys=[source_link_id])

    __tablename__ = 'feature'
    __table_args__ = (
        CheckConstraint(position >= 0),
        CheckConstraint(size > 0),
        {}
    )

    @hybrid_property
    def sequence(self):
        return func.substr(self.component.parent.sequence, self.start, self.end)

    # @cache.memoize()
    def translate(self, component):
        pass

    def __repr__(self):
        return "BoundFeature({})".format(self.id)
