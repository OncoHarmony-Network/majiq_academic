from waitress import serve

from voila.flask.index import app

if __name__ == '__main__':
    serve(app, listen='127.0.0.1:55555')
