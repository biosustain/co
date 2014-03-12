from functools import wraps
from flask import current_app, request
from flask.ext.presst import Relationship, PresstResource, resource_method, fields, ModelResource
from flask.ext.restful import marshal_with, unpack
#from colib.server.models import cache
from colib.mutations import Mutation
from colib.server.models import Component, FeatureProxy, BoundMutation


class marshal_field(object):
    def __init__(self, field):
        """
        :param field: a single field with which to marshal the output.
        """
        self.field = field

    def __call__(self, f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            resp = f(*args, **kwargs)
            if isinstance(resp, tuple):
                data, code, headers = unpack(resp)
                return self.field.format(data), code, headers
            return self.field.format(resp)
        return wrapper


class ComponentResource(ModelResource):

    features = Relationship('FeatureResource')
    mutations = Relationship('MutationResource')

    class Meta:
        model = Component
        exclude_fields = ['sequence']

    def create_item(self, obj):
        pass

    # @resource_method('GET')
    # def features_test(self, component):
    #     features = component.get_features()
    #     return features.all()

    @resource_method('GET')
    @marshal_field(fields.ToMany('feature', embedded=True))
    def added_features(self, component):
        # TODO pagination (important!)
        return map(FeatureProxy.from_record, component.added_features)

    @resource_method('GET')
    @marshal_field(fields.ToMany('feature', embedded=True))
    def removed_features(self, component):
        # TODO pagination (important!)
        return map(FeatureProxy.from_record, component.removed_features)

    @resource_method('POST')
    @marshal_field(fields.ToOne(Component, embedded=True))
    def mutate(self, component, mutations, strict=True):
        """
        :returns: a new component
        """
        new_component = component.mutate(mutations, strict=strict)

        session = self._get_session()
        session.add(new_component)
        session.commit()

        return new_component

    mutate.add_argument('mutations', type=lambda mutations: list(Mutation(**m) for m in mutations))
    mutate.add_argument('strict', type=fields.Boolean(), default=True)

    @resource_method('GET')
    def GET_diff(self, component, other, algorithm=None):
        pass

    @resource_method('GET')
    def GET_fdiff(self, component, other):
        pass

    @resource_method('GET')
    @marshal_field(fields.ToMany(Component))
    def children(self, component):
        return component.descendants

    @resource_method('GET')
    # @cache.memoize()
    @marshal_field(fields.ToMany(Component))
    def lineage(self, component, inclusive):
        return component.get_lineage(inclusive)

    lineage.add_argument('inclusive', location='args', type=lambda b: b.lower() == 'true', default=True)

    @resource_method('GET')
    #@cache.memoize()
    def sequence(self, component):
        return str(component.sequence)

    @resource_method('GET')
    #@cache.memoize()
    def GET_size(self, component):
        return component.size


class FeatureResource(PresstResource):
    type = fields.String()
    name = fields.String()
    position = fields.Integer()
    size = fields.Integer()

    _parent_resource = ComponentResource
    # _parent_id_field = 'component_id'

    class Meta:
        id_field = 'identifier'
        resource_name = 'feature'

    @classmethod
    def get_item_list_for_relationship(cls, relationship, component):
        assert relationship == 'features'
        return component.get_features()

    @resource_method('GET', collection=True)
    def GET_find(self, features, type=None, name=None, qualifiers=None):
        pass

    @resource_method('GET')
    def GET_translate(self, feature, component):
        pass # TODO need to build a chained translation table.

    @resource_method('GET')
    #@cache.memoize()
    def GET_sequence(self, feature):
        return str(feature.sequence)

    @classmethod
    def item_get_resource_uri(cls, item):
        if cls.api is None:
            raise RuntimeError("{} has not been registered as an API endpoint.".format(cls.__name__))
        return '{0}/{1}/{2}'.format(cls.api.prefix, cls.resource_name, item.identifier) # FIXME handle both item and attr.


class MutationResource(ModelResource):

    class Meta:
        resource_name = 'mutation'
        model = BoundMutation


class StrainResource(PresstResource):

    @resource_method('GET')
    @marshal_with(fields.ToMany('StrainResource'))
    def GET_lineage(self, strain):
        return strain.get_lineage()

    @resource_method('GET')
    #@cache.memoize()
    def GET_diff(self, strain, other):
        return strain.diff(other)

    GET_diff.add_argument('other', type=fields.ToOne('StrainResource'))

    @resource_method('GET')
    def GET_tree(self, strains, restrict=None):
        pass


class StrainComponentResource(PresstResource):
    pass


#
# ComponentSchema = Schema.Dict(
#    parent=Schema.Resource('Component', required=False),
#    mutations=Schema.List(Schema.Resource('Mutation'), required=False)
# )

