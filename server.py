from flask import Flask
from flask.ext.presst import PresstApi
from colib.server.models import db
from colib.server.resources import ComponentResource

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'

api = PresstApi(app)
api.add_resource(ComponentResource)

db.init_app(app)
#
# with app.app_context():
#     db.create_all()

if __name__ == "__main__":
    app.run()
