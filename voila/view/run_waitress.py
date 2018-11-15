from waitress import serve

from voila.view.views import app

if __name__ == '__main__':
    serve(app, listen='127.0.0.1:55555')
