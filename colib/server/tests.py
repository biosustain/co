from __future__ import unicode_literals
import logging
import unittest
from flask import json
from flask.ext.testing import TestCase
from colib.mutations import DEL
from colib.server.models import db, Component, BoundFeature
from colib.server.resources import ComponentResource
import server

# logging.basicConfig()
# logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

class ServerTestCase(TestCase):

    def create_app(self):
        server.app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://testuser:testpass@localhost/colibtestdb3'
        # server.app.config['SQLALCHEMY_ECHO'] = True
        return server.app

    def setUp(self):
        with server.app.app_context():
            db.create_all()

            component = Component('AGATATATATACGAGAGCCC')
            component_2 = Component('AGATATATAGAGCCC', parent=component)
            feature_1 = BoundFeature(component=component, type='TATA_box', position=3, size=8)

            db.session.add(component)
            db.session.add(feature_1)
            db.session.add(BoundFeature(component=component, type=None, position=5, size=2))
            #db.session.add(component.mutate([DEL(5, 2)]))
            db.session.add(component_2)

            component_2.removed_features.append(feature_1)

            db.session.add(BoundFeature(component=component_2, type='gene', name='tesT', position=10, size=4))
            db.session.add(Component('AGATATATAGAGCCCCCCC', parent=component_2))

            db.session.commit()

            outer_component_1 = Component('ATGC')
            outer_component_1_feature = BoundFeature(component=outer_component_1, name='DNA-letters', position=1, size=4)
            outer_component_2 = Component('AAAATTTTGGGGCCCC', parent=outer_component_1)

            db.session.add(outer_component_1)
            db.session.add(outer_component_1_feature)
            db.session.add(outer_component_2)

            outer_component_2.removed_features.append(outer_component_1_feature)

            db.session.commit()

    def tearDown(self):
        db.session.remove()
        db.drop_all()

    def test_sequence(self):
        self.assertEqual('AGATATATATACGAGAGCCC', self.client.get('/component/1/sequence').json)
        self.assertEqual('AAAATTTTGGGGCCCC', self.client.get('/component/5/sequence').json)

    def test_size(self):
        self.assertEqual(20, self.client.get('/component/1/size').json)
        self.assertEqual(16, self.client.get('/component/5/size').json)

    def test_children(self):
         self.assertEqual(['/component/2'], self.client.get('/component/1/children').json)
         self.assertEqual([], self.client.get('/component/3/children').json)

    def test_lineage(self):
        self.assertEqual(['/component/1', '/component/2', '/component/3'], self.client.get('/component/3/lineage').json)
        self.assertEqual(['/component/1', '/component/2'], self.client.get('/component/3/lineage?inclusive=False').json)
        self.assertEqual([], self.client.get('/component/1/lineage?inclusive=False').json)

    def test_inherited_features(self):
        self.assertEqual([{'position': 5,
                           'resource_uri': '/feature/1.2t3',
                           'type': None,
                           'name': None,
                           'size': 2},
                          {'position': 10,
                           'resource_uri': '/feature/2.3t3',
                           'type': 'gene',
                           'name': 'tesT',
                           'size': 4}], self.client.get('/component/3/features').json)

    def test_added_features(self):
        self.assertEqual([{'size': 8,
                           'name': None,
                           'resource_uri': '/feature/1.1',
                           'position': 3,
                           'type': 'TATA_box'},
                          {'size': 2,
                           'name': None,
                           'resource_uri': '/feature/1.2',
                           'position': 5,
                           'type': None}],
                         self.client.get('/component/1/added_features').json)

        self.assertEqual([], self.client.get('/component/3/added_features').json)

    def test_mutate(self):
        self.assertEqual('/component/6', self.client.post('/component/1/mutate', data=json.dumps(dict(
            mutations=[DEL(2, 5).__dict__]
        )), content_type='application/json').json['resource_uri'])

        # locked because has children (should have children before):
        self.assertEqual(True, self.client.get('/component/1').json['is_locked'])

        # AGATATATATACGAGAGCCC
        # AG-----TATACGAGAGCCC c.2-5_del
        #    ||||||||          TATA_box  1.2
        #        ||||                    6.5 (or 6.1)

        self.assertEqual('AGTATACGAGAGCCC', self.client.get('/component/6/sequence').json)
        self.assertEqual([{'position': 2,
                           'new_sequence': '',
                           'resource_uri': '/mutation/1',
                           'size': 5}], self.client.get('/component/6/mutations').json)

        self.assertEqual(['/feature/1.1', '/feature/1.2'], self.client.get('/component/6/removed_features').json)

        self.assertEqual([{'size': 10,
                           'position': 2,
                           'resource_uri': '/feature/6.5',
                           'type': 'TATA_box',
                           'name': None}], self.client.get('/component/6/added_features').json)
