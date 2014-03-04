import logging
from flask.ext.testing import TestCase
import sqlalchemy
from colib.mutations import DEL
from colib.server.models import db, Component
import server

logging.basicConfig()
logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

class ServerTestCase(TestCase):

    def create_app(self):
        self.TEST_DB_NAME = 'colibtestdb3'
        server.app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://testuser:testpass@localhost/{}'.format(self.TEST_DB_NAME)
        # server.app.config['SQLALCHEMY_ECHO'] = True
        return server.app

    def setUp(self):
        self.engine = sqlalchemy.create_engine('postgresql://testuser:testpass@localhost/postgres')
        self.connection = self.engine.connect()
        self.connection.execute("COMMIT")
        # self.connection.execute("""SELECT
        #         pg_terminate_backend(procpid)
        #     FROM
        #         pg_stat_activity
        #     WHERE procpid <> pg_backend_pid() AND datname = '{}'""".format(self.TEST_DB_NAME))
        self.connection.execute("DROP DATABASE IF EXISTS {}".format(self.TEST_DB_NAME))
        self.connection.execute("COMMIT")
        self.connection.execute("CREATE DATABASE {}".format(self.TEST_DB_NAME))

        with server.app.app_context():
            db.create_all()

            component = Component('AGATATACGAGAGCCC')
            component_2 = Component('AGATAGAGCCC', parent=component)


            db.session.add(component)
            #db.session.add(component.mutate([DEL(5, 2)]))
            db.session.add(component_2)
            db.session.add(Component('AGATAGAGCCCCCCC', parent=component_2))
            db.session.commit()

            print(Component.query.all())

        print(server.app.url_map)


    def test_lineage(self):
        self.assertEqual(['/component/1', '/component/2', '/component/3'], self.client.get('/component/3/lineage').json)
        self.assertEqual(['/component/1', '/component/2'], self.client.get('/component/3/lineage?inclusive=False').json)
        self.assertEqual([], self.client.get('/component/1/lineage?inclusive=False').json)

    #
    # def test_inherited_features(self):
    #     self.assertEqual(None, self.client.get('/component/2/features_test').json)
    #
    # def tearDown(self):
    #     self.connection.execute("COMMIT")
    #     self.connection.execute("DROP DATABASE {}".format(self.TEST_DB_NAME))