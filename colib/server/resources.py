from functools import wraps
from flask import current_app
from flask.ext.presst import Relationship, PresstResource, resource_method, fields, ModelResource
from flask.ext.restful import marshal_with, unpack
#from colib.server.models import cache
from colib.server.models import Component

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

    #features = Relationship('FeatureResource')
    mutations = Relationship('MutationResource')

    class Meta:
        model = Component

    def create_item(self, obj):
        pass

    @resource_method('GET')
    def features_test(self, component):
        features = component.get_features()
        print(features.all())
        return []

    @resource_method('GET')
    @marshal_with(fields.ToMany('FeatureResource'))
    def GET_added_features(self, component):
        return component.features

    @resource_method('GET')
    @marshal_field(fields.ToMany('FeatureResource'))
    def GET_removed_features(self, component):
        return component.removed_features

    @resource_method('POST')
    @marshal_field(fields.ToOne('ComponentResource'))
    def POST_mutate(self, component, mutations, strict=True):
        """
        :returns: a new component
        """
        new_component = component.mutate(mutations, strict=strict)

        return new_component

    POST_mutate.add_argument('mutations', type=fields.ToMany('MutationResource'))
    POST_mutate.add_argument('strict', type=fields.Boolean(), default=True)

    @resource_method('GET')
    def GET_diff(self, component, other, algorithm=None):
        pass

    @resource_method('GET')
    def GET_fdiff(self, component, other):
        pass

    @resource_method('GET')
    #@cache.memoize()
    #@marshal_field(fields.ToMany('ComponentResource'))
    def lineage(self, component, inclusive):
        print(component.get_lineage())
        return [self.item_get_resource_uri(c) for c in component.get_lineage(inclusive)]

    lineage.add_argument('inclusive', location='args', type=lambda b: b == 'True', default=True)

    @resource_method('GET')
    #@cache.memoize()
    def GET_sequence(self, component):
        return str(component.sequence)

    @resource_method('GET')
    #@cache.memoize()
    def GET_size(self, component):
        return component.size


class FeatureResource(PresstResource):

    class Meta:
        resource_name = 'feature'

    def get_item_list_for_relationship(cls, relationship, parent_item):
        # TODO check request query for include_inherited.

        pass

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


class MutationResource(PresstResource):
    pass


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

