import sys, os
os.chdir(os.path.dirname(__file__))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/python_libs')

import inzyme
application = inzyme.wsgi_app()

