import os.path

import betastasis_data as bdata

# Fall back to StringIO in environments where cStringIO is not available
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from django.conf import settings
from django.template.loader import get_template
from django.template import Context
from django.template import RequestContext
from django.http import HttpResponse
from django.http import Http404
from django.shortcuts import render_to_response
from django.utils import simplejson

def app_disabled_view(request):
    return HttpResponse('Nothing here')

def cancer_type_list(request):
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    return render_to_response('cancer_types.html',
                              {'cancer_data_sets': cancer_data_sets},
                              context_instance=RequestContext(request))

def cancer_type_data_sets(request, cancer_type):
    if cancer_type not in bdata.cancer_types:
        raise Http404
    
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    context = {'current_cancer': bdata.make_current_cancer_dict(cancer_type),
               'cancer_data_sets': cancer_data_sets}
    
    return render_to_response('data_set_list.html',
                              context,
                              context_instance=RequestContext(request))

def set_platform_list(request, cancer_type, data_set_id):
    if data_set_id not in bdata.data_sources:
        raise Http404
    
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    context = {'current_cancer': bdata.make_current_cancer_dict(cancer_type),
               'current_data_set': {'id': data_set_id, 'meta': bdata.data_sources[data_set_id]},
               'cancer_data_sets': cancer_data_sets}
    
    return render_to_response('data_set_platforms.html',
                              context,
                              context_instance=RequestContext(request))
    
def platform_gene_list_json(request, data_set_id, platform_id):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'expression')
    if path is None:
        raise Http404
    
    # Compose the path to the gene list .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, path)
    gene_list_path = os.path.join(platform_root, "genes.json")
    if not os.path.isfile(gene_list_path):
        raise Http404
    
    f = open(gene_list_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')
    
def expression_json(request, data_set_id, platform_id, gene):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'expression')
    if path is None:
        raise Http404
    
    # Compose the path to the expression data .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, path)
    gene_path = os.path.join(platform_root, gene[0].lower())
    gene_path = os.path.join(gene_path, gene.upper() + ".json")
    if not os.path.isfile(gene_path):
        raise Http404
    
    f = open(gene_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')

def expression_boxplot(request, cancer_type, data_set_id, platform_id):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'expression')
    if path is None:
        raise Http404
    
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    context = {'current_cancer': bdata.make_current_cancer_dict(cancer_type),
               'current_data_set': {'id': data_set_id, 'meta': bdata.data_sources[data_set_id]},
               'cancer_data_sets': cancer_data_sets,
               'current_platform': {'id': platform_id, 'name': bdata.platforms[platform_id]},
               'current_data_type': 'expression',
               'default_gene': 'TP53'}
    
    return render_to_response('gene_boxplot.html',
                              context,
                              context_instance=RequestContext(request))

def platform_survival_gene_list_json(request, data_set_id, platform_id):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'survival')
    if path is None:
        raise Http404

    # Compose the path to the gene list .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, path)
    gene_list_path = os.path.join(platform_root, "genes.json")
    if not os.path.isfile(gene_list_path):
        raise Http404
    
    f = open(gene_list_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')
    
def survival_expression_json(request, data_set_id, platform_id, gene):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'survival')
    if path is None:
        raise Http404
    
    # Compose the path to the expression data .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, path)
    gene_path = os.path.join(platform_root, gene[0].lower())
    gene_path = os.path.join(gene_path, gene.upper() + ".json")
    if not os.path.isfile(gene_path):
        raise Http404
    
    f = open(gene_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')

def survival_json(request, data_set_id, platform_id):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'survival')
    if path is None:
        raise Http404
    
    # Compose the path to the survival data .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, path)
    gene_path = os.path.join(platform_root, "survival.json")
    if not os.path.isfile(gene_path):
        raise Http404
    
    f = open(gene_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')

def kaplan_meier_plot(request, cancer_type, data_set_id, platform_id):
    path = bdata.get_data_path_if_exists(data_set_id, platform_id, 'survival')
    if path is None:
        raise Http404
    
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    
    context = {'current_cancer': bdata.make_current_cancer_dict(cancer_type),
               'current_data_set': {'id': data_set_id, 'meta': bdata.data_sources[data_set_id]},
               'cancer_data_sets': cancer_data_sets,
               'current_platform': {'id': platform_id, 'name': bdata.platforms[platform_id]},
               'current_data_type': 'survival'}
    
    return render_to_response('gene_survival.html',
                              context,
                              context_instance=RequestContext(request))

def features_json(request, cancer_type, data_set_id, platform_id):
    file_name = bdata.get_data_path_if_exists(data_set_id, platform_id, 'feature_matrix')
    if file_name is None:
        raise Http404
    
    # Compose the path to the .json-file
    features_path = os.path.join(settings.DATA_ROOT_PATH, file_name)
    if not os.path.isfile(features_path):
        raise Http404
    
    f = open(features_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')    

def feature_matrix_plot(request, cancer_type, data_set_id, platform_id):
    file_name = bdata.get_data_path_if_exists(data_set_id, platform_id, 'feature_matrix')
    if file_name is None:
        raise Http404
    
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    
    context = {'current_cancer': bdata.make_current_cancer_dict(cancer_type),
               'current_data_set': {'id': data_set_id, 'meta': bdata.data_sources[data_set_id]},
               'cancer_data_sets': cancer_data_sets,
               'current_platform': {'id': platform_id, 'name': bdata.platforms[platform_id]},
               'current_data_type': 'feature_matrix'}
    
    return render_to_response('feature_matrix.html',
                              context,
                              context_instance=RequestContext(request))

def probeset_genelist_json(request, platform_id, transcriptome):
    path = bdata.get_probeset_path_if_exists(platform_id, transcriptome)
    print path
    if path is None:
        raise Http404

    # Compose the path to the gene list .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, 'platforms')
    platform_root = os.path.join(platform_root, path)
    gene_list_path = os.path.join(platform_root, "genelist.json")
    if not os.path.isfile(gene_list_path):
        raise Http404
    
    f = open(gene_list_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')

def probeset_gene_json(request, platform_id, transcriptome, gene):
    path = bdata.get_probeset_path_if_exists(platform_id, transcriptome)
    if path is None:
        raise Http404

    # Compose the path to the .json-file
    platform_root = os.path.join(settings.DATA_ROOT_PATH, 'platforms')
    platform_root = os.path.join(platform_root, path)
    gene_path = os.path.join(platform_root, gene[0].lower())
    gene_path = os.path.join(gene_path, gene.upper() + ".json")
    if not os.path.isfile(gene_path):
        raise Http404
    
    f = open(gene_path, 'rb')
    
    return HttpResponse(f, mimetype='application/json')

def probeset_design_plot(request, platform_id, transcriptome):
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    probesets = bdata.get_probesets_by_platform()
    
    context = {'current_platform': {'id': platform_id, 'name': bdata.platforms[platform_id]},
               'cancer_data_sets': cancer_data_sets,
               'transcriptome': {'id': transcriptome, 'name': bdata.transcriptomes[transcriptome]},
               'probesets': probesets,
               'current_data_type': 'probeset_design'}
    
    return render_to_response('probeset_design.html',
                              context,
                              context_instance=RequestContext(request))

def all_platforms_list(request):
    cancer_data_sets = bdata.get_data_sets_by_cancer_type()
    probesets = bdata.get_probesets_by_platform()
    
    context = {'cancer_data_sets': cancer_data_sets,
               'probesets': probesets}
    
    return render_to_response('all_platforms.html',
                              context,
                              context_instance=RequestContext(request))
