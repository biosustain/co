from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy import func, CheckConstraint, select, not_
from sqlalchemy.dialects import postgres
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.mutable import MutableDict
from sqlalchemy.orm import backref, deferred, aliased
from werkzeug.utils import cached_property
import colib

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
    removed_features = db.relationship('BoundFeature',
                                       secondary=component_removed_features_map,
                                       backref='removed_from_components')

    is_locked = db.Column(db.Boolean(), default=False)

    def __init__(self, sequence=None, parent=None):
        self.sequence = sequence
        self.parent = parent

    @hybrid_property
    def is_mutated(self):
        return self.parent is not None

    @hybrid_property
    def size(self):
        return len(self.sequence)

    @size.expression
    def size(self):
        return func.length(self.sequence)

    def get_features(self, include_inherited=True):
        """
        This function executes an SQL query roughly equivalent to:

            WITH RECURSIVE parent_components(id, parent_id) AS
            (SELECT component.id AS id, component.parent_id AS parent_id
            FROM component
            WHERE component.id = 3 UNION ALL SELECT c.id AS c_id, c.parent_id AS c_parent_id
            FROM component AS c, parent_components AS a
            WHERE c.id = a.parent_id)
             SELECT *
            FROM feature
            INNER JOIN parent_components ON (parent_components.id = feature.component_id)
            WHERE feature.component_id = parent_components.id
                AND NOT EXISTS (SELECT 1 FROM component_removed_features_map WHERE
                    component_removed_features_map.component_id IN (SELECT id FROM parent_components) AND
                    component_removed_features_map.feature_component_id = feature.component_id AND
                    component_removed_features_map.feature_id = feature.id)
            ORDER BY feature.position
        """
        if include_inherited:
            parent_components = db.session.query(Component.id, Component.parent_id)\
                .filter(Component.id==self.id)\
                .cte(name="parent_components", recursive=True)

            parents_alias = aliased(parent_components, name='a')
            components_alias = aliased(Component, name="c")
            parent_components = parent_components.union_all(
                db.session.query(components_alias.id, components_alias.parent_id)\
                    .filter(components_alias.id==parents_alias.c.parent_id))

            deleted_from_ancestor = db.session.query(component_removed_features_map)\
                .filter(component_removed_features_map.c.component_id.in_(select([parent_components.c.id])),
                        component_removed_features_map.c.feature_component_id == BoundFeature.component_id,
                        component_removed_features_map.c.feature_id == BoundFeature.id).exists()

            query = BoundFeature.query\
                .filter(BoundFeature.component_id == parent_components.c.id)\
                .filter(not_(deleted_from_ancestor))\
                .order_by(BoundFeature.position)

            # TODO need to return some kind of feature-queryset proxy that converts to a FeatureProxy only when the result is loaded.
            return (FeatureProxy.from_record(feature, self) for feature in query)
        else:
            return self.added_features
        # TODO map.

    def get_lineage(self, inclusive=True):
        """

        """
        ancestors = db.session.query(Component.id, Component.parent_id)\
            .filter(Component.id==self.id)\
            .cte(name="ancestors", recursive=True)

        ancestors_alias = aliased(ancestors, name='a')
        components_alias = aliased(Component, name="c")
        ancestors = ancestors.union_all(
            db.session.query(components_alias.id, components_alias.parent_id)
                .filter(components_alias.id == ancestors_alias.c.parent_id))

        filter_id = ancestors.c.id if inclusive else ancestors.c.parent_id
        return Component.query.filter(Component.id == filter_id)

    def find_features(self, type=None, name=None, qualifiers=None, xref=None):
        pass

    def _make_co_object(self, include_parent=False):

        if include_parent and self.parent is not None:
            parent = self.parent.get_component_obj()
        else:
            parent = None

        co = colib.Component(sequence=self.sequence, parent=parent, feature_class=FeatureProxy)

        # TODO add features attribute to colib.Component()
        for feature in self.get_features(include_inherited=not include_parent):
            co.features.add(feature)

        return co

    @cached_property
    def co_object(self):
        return self._make_co_object()

    def mutate(self, mutations, strict=True):
        """

        """
        colib_mutated = self.co_object.mutate(mutations, strict=strict)

        component = Component(sequence=str(colib_mutated), parent=self)

        db.session.add(component)

        for feature in colib_mutated.features.removed:
            component.removed_features.append(feature.record)

        for feature in colib_mutated.features.added:
            component.added_features.append(feature.record.move(feature.position, feature.size, component))

        for mutation in mutations:
            component.mutations.append(BoundMutation(**mutation.__dict__))

        self.is_locked = True
        return component

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


class FeatureProxy(colib.Feature):
    def __init__(self, component, position, size, record=None, record_id=None, strand=None, type=None, name=None, translated_view=None):
        super(FeatureProxy, self).__init__(component, position, size, strand=None, type=type, name=name)
        self.record = record
        self.record_id = record_id
        self.translated_view = translated_view

    @classmethod
    def from_record(cls, feature, component=None):
        if component is not None and feature.component != component:
            return FeatureProxy(
                feature.component,
                feature.position,
                feature.size,
                feature,
                feature.id,
                type=feature.type,
                name=feature.name,
                translated_view=component)

        return FeatureProxy(
            feature.component,
            feature.position,
            feature.size,
            feature,
            feature.id,
            type=feature.type,
            name=feature.name)

    @property
    def sequence(self):
        sequence = self._component.co_object[self.start:self.end + 1]
        if self.strand == colib.REVERSE_STRAND:
            return sequence.reverse_complement()
        return sequence

    def translate_to(self, component, using_tt=None):
        if component == self.component.co_object:
            return self
        if using_tt is None:
            using_tt = component.tt(self.component.co_object)
        print(self._position, using_tt[self.position])
        return self.move(using_tt[self.position], size=self.size, component=component)

    @property
    def identifier(self):
        if self.translated_view:
            return '{}.{}t{}'.format(self.component.id, self.record_id, self.translated_view.id)
        else:
            return '{}.{}'.format(self.component.id, self.record_id)

    def __repr__(self):
        return 'FeatureProxy({}: {} to {})'.format(self.type, self.start, self.end)


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

    strand = db.Column(db.Enum('forward', 'reverse', name='strand'))

    type = db.Column(db.String, index=True)  # Sequence Ontology terms strongly preferred.
    name = db.Column(db.String, index=True)
    #
    #meta = db.Column(postgres.JSON, nullable=False)

    source_id = db.Column(db.Integer, db.ForeignKey(Component.id))
    source_sequence_is_match = db.Column(db.Boolean, default=True)

    component = db.relationship(Component,
                                foreign_keys=[component_id],
                                backref=backref('added_features', lazy='dynamic', order_by=position))

    source = db.relationship(Component, foreign_keys=[source_id])

    __tablename__ = 'feature'
    __table_args__ = (
        CheckConstraint(position >= 0),
        CheckConstraint(size > 0),
        {}
    )

    def move(self, position, size, component=None):
        """
        Returns a copy of the feature with the new position, size, and component.
        """
        return BoundFeature(
            component=component,
            position=position,
            size=size,
            strand=self.strand,
            type=self.type,
            source=self.source,
            source_sequence_is_match=False # TODO compute source sequence match.
        )

    @hybrid_property
    def sequence(self):
        return func.substr(self.component.parent.sequence, self.start, self.end)

    # @cache.memoize()
    def translate(self, component):
        pass

    def __repr__(self):
        return "BoundFeature({})".format(self.id)

