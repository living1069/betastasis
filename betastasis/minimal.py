import os, sys
#sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/python_libs')

import bottle
from bottle import route

bottle.debug(True)

def wsgi_app():
	return bottle.default_app()

@route('/')
def index():
	return '<h1>Hello world!</h1>'

