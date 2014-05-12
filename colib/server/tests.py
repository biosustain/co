from __future__ import unicode_literals

from flask import json
from flask.ext.testing import TestCase

from colib.mutations import DEL, INS
from colib.server.models import db, Component, BoundFeature
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
            component.added_features.append(BoundFeature(type='TATA_box', position=3, size=8))
            component.added_features.append(BoundFeature(type=None, position=5, size=2))
            db.session.add(component)
            db.session.commit()

            component_2 = component.mutate([DEL(5, 4)])
            component_2.added_features.append(BoundFeature(type='gene', name='tesT', position=10, size=4))
            db.session.add(component_2)
            db.session.commit()

            # component_2.removed_features.append(feature_1)

            db.session.add(component_2.mutate([INS(15, 'CCCC')]))
            db.session.commit()

            outer_component_1 = Component('ATGC')
            outer_component_1_feature = BoundFeature(component=outer_component_1, name='DNA-letters', position=1,
                                                     size=4)
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
        self.assertEqual('AGATATACGAGAGCCC', self.client.get('/component/2/sequence').json)
        self.assertEqual('AGATATACGAGAGCCCCCCC', self.client.get('/component/3/sequence').json)
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
        self.assertEqual([{'type': 'TATA_box',
                           'name': None,
                           'resource_uri': '/feature/1.1',
                           'size': 8,
                           'position': 3},
                          {'name': None,
                           'position': 5,
                           'resource_uri': '/feature/1.2',
                           'size': 2,
                           'type': None}], self.client.get('/component/1/features').json)
        # What is feature 3?

        print(self.client.get('/component/2/features').json)
        self.assertEqual([{'type': 'TATA_box',
                           'size': 4,  # mid section deleted
                           'resource_uri': '/feature/2.3',
                           'position': 3,
                           'name': None},
                          {'type': 'gene',
                           'size': 4,
                           'resource_uri': '/feature/2.4',
                           'position': 10,
                           'name': 'tesT'}], self.client.get('/component/2/features').json)

        self.assertEqual([{'name': None,
                           'position': 3,
                           'resource_uri': '/feature/2.3t3',
                           'size': 4,
                           'type': 'TATA_box'},
                          {'type': 'gene',
                           'size': 4,
                           'resource_uri': '/feature/2.4t3',
                           'position': 10,
                           'name': 'tesT'}], self.client.get('/component/3/features').json)

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
                           'resource_uri': '/mutation/3',
                           'size': 5}], self.client.get('/component/6/mutations').json)

        self.assertEqual(['/feature/1.1', '/feature/1.2'], self.client.get('/component/6/removed_features').json)

        self.assertEqual([{'size': 4,
                           'position': 2,
                           'resource_uri': '/feature/6.6',
                           'type': 'TATA_box',
                           'name': None}], self.client.get('/component/6/added_features').json)
