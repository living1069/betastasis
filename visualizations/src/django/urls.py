from django.conf.urls.defaults import *
from django.conf import settings

from betastasis.views import *

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

betastasis_patterns = patterns('',
	('^$', cancer_type_list),
	('^platforms/$', all_platforms_list),
	('^platforms/(?P<platform_id>\w+)/(?P<transcriptome>\w+)/gene_list_json/$', probeset_genelist_json),
	('^platforms/(?P<platform_id>\w+)/(?P<transcriptome>\w+)/(?P<gene>\w+)/gene_json/$', probeset_gene_json),
	('^platforms/(?P<platform_id>\w+)/(?P<transcriptome>\w+)/probeset_design/$', probeset_design_plot),
	('^(?P<cancer_type>\w+)/$', cancer_type_data_sets),  
	('^(?P<cancer_type>\w+)/(?P<data_set_id>\w+)/$', set_platform_list),
	('^(?P<cancer_type>\w+)/(?P<data_set_id>\w+)/(?P<platform_id>\w+)/expr_boxplot/$', expression_boxplot),
	('^(?P<cancer_type>\w+)/(?P<data_set_id>\w+)/(?P<platform_id>\w+)/survival_plot/$', kaplan_meier_plot),
	('^(?P<cancer_type>\w+)/(?P<data_set_id>\w+)/(?P<platform_id>\w+)/feature_matrix/$', feature_matrix_plot),
	('^(?P<data_set_id>\w+)/(?P<platform_id>\w+)/gene_list_json/$', platform_gene_list_json),
	('^(?P<data_set_id>\w+)/(?P<platform_id>\w+)/(?P<gene>\w+)/expr_json/$', expression_json),
	('^(?P<data_set_id>\w+)/(?P<platform_id>\w+)/surv_json/$', survival_json),
	('^(?P<data_set_id>\w+)/(?P<platform_id>\w+)/(?P<gene>\w+)/surv_expr_json/$', survival_expression_json),
	('^(?P<cancer_type>\w+)/(?P<data_set_id>\w+)/(?P<platform_id>\w+)/features_json/$', features_json),
)

urlpatterns = None

if settings.APP_TOGGLE is False:
	urlpatterns = patterns('',
			       ('^$', app_disabled_view))
else:
	urlpatterns = betastasis_patterns

# Set up 
if settings.DEBUG:
	media_url = settings.MEDIA_URL[1:]
	static_pattern = r'^%s(?P<path>.*)$' % media_url
	urlpatterns += patterns('',
				(static_pattern, 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT,
									       'show_indexes': True}),
                                )
