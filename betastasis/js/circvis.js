function drawCircvis (data) {
        parseData(data);
        wedge_plot(document.getElementById('fig'), parsed_data);
};

function parseData(data) {
    parsed_data = { complete_network: [], features:[], located_features:[], network:[], unlocated:[], unlocated_features:[]};
    //console.log(data);

    var data_array = data.split('\n');
    var node_index = 0;
    for (var i=0; i < data_array.length-1; i++)
    {
        var row_array = data_array[i].split('\t');
        var feature1 = row_array[0].split(':');
        var feature2 = row_array[1].split(':');
        var color = null;
        if (row_array.length > 2){
            color = row_array[2];
        }


        var end = parseInt(feature1[2]) || 0;

        var node1={label:row_array[0], chr:feature1[0], start:feature1[1], end:feature1[2], index:node_index, level:0, nodeName:':'+row_array[0].replace(/:/g,''), source:"", tile_length:end - parseInt(feature1[1]), linkDegree:0};

        ++node_index;
        var end = parseInt(feature2[2]) || 0;
        var node2={label:row_array[1], chr:feature2[0], start:feature2[1], end:feature2[2], index:node_index, level:0, nodeName:':'+row_array[1].replace(/:/g,''), source:"", tile_length:end - parseInt(feature2[1]), linkDegree:0};
        ++node_index;
        
        parsed_data.complete_network.push({'node1':node1,'node2':node2, value:color, yscore:color});
        parsed_data.network.push({'node1':node1,'node2':node2, value:color, yscore:color});
    }
};


function wedge_plot(div, parsed_data) {
    var width=600, height=600;
    var	ring_radius = width / 20;
    
    var nodeColor = {"GEXP":"red","GENO":"pink","PROT":"cornflowerblue","METAB":"white","PHENO":"lightgrey", "None":"black"};

    var chrom_keys = ["1","2","3","4","5","6","7","8","9","10",
        "11","12","13","14","15","16", "17","18","19","20","21","22","X","Y"];

    var annotations= {chrom_leng: [{'chr_name':'1','chr_length':249250621},
					{'chr_name':'2','chr_length':243199373},
					{'chr_name':'3','chr_length':198022430},
					{'chr_name':'4','chr_length':191154276},
					{'chr_name':'5','chr_length':180915260},
					{'chr_name':'6','chr_length':171115067},
					{'chr_name':'7','chr_length':159138663},
					{'chr_name':'8','chr_length':146364022},
					{'chr_name':'9','chr_length':141213431},
					{'chr_name':'10','chr_length':135534747},
					{'chr_name':'11','chr_length':135006516},
					{'chr_name':'12','chr_length':133851895},
					{'chr_name':'13','chr_length':115169878},
					{'chr_name':'14','chr_length':107349540},
					{'chr_name':'15','chr_length':102531392},
					{'chr_name':'16','chr_length':90354753},
					{'chr_name':'17','chr_length':81195210},
					{'chr_name':'18','chr_length':78077248},
					{'chr_name':'19','chr_length':59128983},
					{'chr_name':'20','chr_length':63025520},
					{'chr_name':'21','chr_length':48129895},
					{'chr_name':'22','chr_length':51304566},
					{'chr_name':'X','chr_length':155270560},
					{'chr_name':'Y','chr_length':59373566}
					]};

	var annotation_ring = [];
	
    //var chrom_translate_hash = {"1":"chr01","2":"chr02","3":"chr03","4":"chr04","5":"chr05","6":"chr06","7":"chr07","8":"chr08","9":"chr09","10":"chr10",
      //  "11":"chr11","12":"chr12","13":"chr13","14":"chr14","15":"chr15","16":"chr16", "17":"chr17", "18":"chr18", 
       // "19":"chr19", "20":"chr20", "21":"chr21", "22":"chr22", "X":"chrX", "Y":"chrY"};	

   /* var chrom_keys = ["I","II","III","IV","V","VI","VII","VIII","IX","X",
        "XI","XII","XIII","XIV","XV","XVI", "MT"];
    var chrom_translate_hash = {"I":"chr01","II":"chr02","III":"chr03","IV":"chr04","V":"chr05","VI":"chr06","VII":"chr07","VIII":"chr08","IX":"chr09","X":"chr10",
        "XI":"chr11","XII":"chr12","XIII":"chr13","XIV":"chr14","XV":"chr15","XVI":"chr16", "MT":"chrmt"};	*/
    function genome_listener(chr) {
            //renderLinearData(chr);
        }
    function wedge_listener(feature) {
                    var chr = feature.chr;
                    var start = bpToMb(feature.start) - 2.5;
                    var range_length = 100000;//bpToMb(feature.end) - start + 2.5;
                    renderLinearData(chr,start,linear_unit);
    }
    function reference_tile_listener(feature) {                    
		    var chr = feature.chr;
		    window.open('http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=chr' + chr + ':' + feature.start + '..' + (parseInt(feature.start) + sgdbrowser_bin),'_blank');
    }
    function node_tile_listener(feature) {
                    var chr = feature.chr;
                    window.open('http://browse.yeastgenome.org/fgb2/gbrowse/scgenome/?name=chr' + chr + ':' + feature.start + '..' + (parseInt(feature.start) + sgdbrowser_bin),'_blank');
    }

    /*function gcn4_listener(feature) {
	var fs = parseInt(feature.start) + 600;
	var fe = parseInt(feature.end) + 600;
	var fchr = chrom_translate_hash[feature.chr];
	window.open(gbrowse_url + 'ref=' + fchr + '&start=' + (fs) +
                      '&end='+ (feature.end ? fe : feature.start + linear_unit),'_blank');
    }*/
     //var ucsc_genome_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=S.+cerevisiae&db=sacCer2&';
          var tick_listener = function(feature){
              window.open(sgd_url + feature.systematic_name,'_blank');
		return false;
              };
       var karyotype_tooltip_items = {
           'Karyotype Label' : function(feature) { return  vq.utils.VisUtils.options_map(feature)['label'];},
            Location :  function(feature) {
				var offset = 8000;
				if (feature.type == null || feature.type == 'uniprobe') 
					offset = 0
				return 'Chr' + feature.chr + ' ' + feature.start + '-' + (parseInt(feature.end) - offset);
			}
        },
	gcn4_tooltip_items = {
           'GCN4' : function(feature) { return  vq.utils.VisUtils.options_map(feature)['label'];},
            Location :  function(feature) {
                                var offset = 8000;
                                if (feature.type == null || feature.type == 'uniprobe')
                                        offset = 0
                                return 'Chr' + feature.chr + ' ' + feature.start + '-' + (parseInt(feature.end) - offset);
                        }
        },
	reference_tooltip_items = {
            "Gene":  function(feature) {
                                //if (feature.label && feature.label.indexOf("hotspot") != -1)
				return feature.label.toUpperCase() + ' ' + feature.chr + ' ' + feature.start;
				//return 'Gene@Chr' + feature.chr + ' ' + feature.start;
                        }
        },
        annotation_tooltip_items = {
            "Genomic Annotation":  function(feature) {
                                return feature.note;
                        }
        },
        unlocated_tooltip_items = {
            "HotRegion" :  function(feature) {
			if (feature.sourceNode.source == 'ALL')
				feature.sourceNode.source = 'Cumulative'; 
			return feature.sourceNode.label + " chr" + feature.sourceNode.chr + ":" + feature.sourceNode.start + "-" + feature.sourceNode.end +" for conditions:" + feature.sourceNode.source;
        }
	};
    var chrom_leng = vq.utils.VisUtils.clone(annotations['chrom_leng']);
    //parsed_data = {network : null,unlocated : null, features : null,unlocated_features:null,located_features:null};
    	var ring = vq.utils.VisUtils.clone(parsed_data['features']).filter(function(feature) { return feature.chr != '' && feature.rank != '';})
           .map(function(feature) { return vq.utils.VisUtils.extend(feature,{value:feature.rank, options:feature.weight});});
    //var weight_color = pv.Scale.linear(ring.map(function(feature) { return feature.options;})).range('blue','red');

    	var ticks = vq.utils.VisUtils.clone(parsed_data['features']);

    var unlocated_map = vq.utils.VisUtils.clone(parsed_data['unlocated']).filter(function(link) { 
	return  link.node1.chr != '';})
            .map(function(link) {
      var node =  vq.utils.VisUtils.extend(link.node2,{ chr:link.node1.chr, start:link.node1.start,end:link.node1.end, value: 0});
        node.sourceNode = vq.utils.VisUtils.extend({},link.node1); node.targetNode = vq.utils.VisUtils.extend({},link.node2);
        node.importance = link.importance;
	return node;
    }).concat(vq.utils.VisUtils.clone(parsed_data['unlocated']).filter(function(link) { return  link.node2.chr != '';})
            .map(function(link) {
      var node =  vq.utils.VisUtils.extend(link.node1,{ chr:link.node2.chr, start:link.node2.start,end:link.node2.end, value: 0});
        node.sourceNode = vq.utils.VisUtils.extend({},link.node1); node.targetNode = vq.utils.VisUtils.extend({},link.node2);
        node.importance = link.importance;
	return node;
    }));
    
var data = {
        GENOME: {
            DATA:{
                key_order : chrom_keys,
                key_length : chrom_leng
            },
            OPTIONS: {
                radial_grid_line_width: 1,
                label_layout_style : 'clock',
                listener : genome_listener,
                label_font_style : '12pt helvetica bold'
            }
        },
        TICKS : {
            DATA : {
                data_array : ticks
            },
            OPTIONS :{
                 display_legend : false,
                listener : wedge_listener,
				overlap_distance : 12000,
                fill_style : function(tick) {
			//return "red";
			return node_colors(tick.source); 
		},
                tooltip_items : {Tick : function(node) { return node.label+ ' ' + node.source + ' Chr' + node.chr + ' ' + node.start +
                            '-' + node.end;}}
            }
        },
        PLOT: {
            width : width,
            height :  height,
            horizontal_padding : 30,
            vertical_padding : 30,
            container : div,
            enable_pan : false,
            enable_zoom : false,
            show_legend: true,
            legend_include_genome : true,
            legend_corner : 'ne',
            legend_radius  : width / 15
        },
           WEDGE:[
          {     
                PLOT : {
                    height : ring_radius/2,
                    type :   'karyotype'
                },
                DATA:{
                    data_array : cytoband
                },
                OPTIONS: {
                    legend_label : '' ,
                    legend_description : 'markers',
                    outer_padding : 2,
                    tooltip_items : reference_tooltip_items,
                    listener : reference_tile_listener
                }
            },  {
                PLOT : {
                    height : ring_radius/2,
                    type :   'karyotype'
                },
                DATA:{
                    data_array : annotation_ring
                },
                OPTIONS: {
                    legend_description : 'Annotation',
                    outer_padding : 5,
                    tooltip_items : annotation_tooltip_items,
                    listener : wedge_listener
                }
            }
        ],
        NETWORK:{
            DATA:{
                data_array : parsed_data['network']
            },
            OPTIONS: {
                outer_padding : 15,
                node_highlight_mode : 'isolate',
                link_line_width : 2,
                node_key : function(node) { 
			var id= node['source'] + ':' + node['chr']  + node["start"];
			return id;
			},
                node_fill_style: function(ynode){
                        if (ynode == undefined || ynode.source == undefined){
                                return "blue";
                        }
                        var col = nodeColor[ynode.source];
                        if (col == undefined){
                                return "blue";
                        }else{
                                return col;
                        }
                },
                node_stroke_style: function() { return 'none'; },
                node_listener : node_tile_listener,
		//node_stroke_style: function(ynode){
	        //	return "blue";
                //},
                //link_listener: initiateScatterplot,
                link_stroke_style : function(link) {
        		   /* var novellink = link.sourceNode.label.split(":")[0] + "_" + link.targetNode.label.split(":")[0];      		
        		    if (novelSet[novellink] == 1 && yanexContext != ""){
        		    	if (yanexContext == "Known"){
        				return "black";
        			}
        			if (yanexContext == "Novel"){
        				return "red";
        			}
        		    }*/
        		    var ycolor = link.value;
        		    if (ycolor == null){
                        return "skyblue";    
        		    }
        		    return ycolor;
        		},
                constant_link_alpha : 0.5,
				tile_nodes : true,
                node_overlap_distance : 1000,
                node_tooltip_items :  {Node : function(node) { return node.label+ ' ' + node.source + ' Chr' + node.chr + ' ' + node.start +
                            '-' + node.end;}},
                link_tooltip_items :  {
                    'Association' : function(link) { return link.sourceNode.label+ ' ' + link.sourceNode.source + ' Chr' + link.sourceNode.chr + ':' + link.sourceNode.start + '-' + link.targetNode.label+ ' ' + link.targetNode.source + ' Chr ' + link.targetNode.chr + ":" + link.targetNode.start; }
                }
            }
        }
    };
    var circle_vis = new vq.CircVis();
    var dataObject ={DATATYPE : "vq.models.CircVisData", CONTENTS : data };
    circle_vis.draw(dataObject);

   if (circle_vis.pan_enable != null) {
       re.circvis_obj.setPanEnabled(circle_vis.pan_enable);
   }
   if (circle_vis.zoom_enable  != null) {
       re.circvis_obj.setZoomEnabled(circle_vis.zoom_enable);
   }    

    return circle_vis;
};