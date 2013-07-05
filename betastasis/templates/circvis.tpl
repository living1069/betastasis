%rebase base path=path, url=url, scripts='ui_circvis,circvis,chrom_data'

%# Here we figure out the dataset title by looking for the "title" attribute.
%_, _, attr = hierarchy['/'.join(url.split('/')[:-1])]
<h2>{{attr['title'] if 'title' in attr else 'Tools'}}</h2>

<h3>{{title if 'title' in locals() else 'Circvis'}}</h3>

<div id="fig" style="width: 600px;"></div>

<div id="sidebar">
    <div class="sidebox">
    %if not data:
         <p>You can upload your own tab separated file to display associations. <br> 
            Example: <br>
            1:59252019::    5:174769::  #000080
            1:62517091::    1:858629:: </p>

		<br>
        <b>Upload</b><br>
            <button id="uploader"><input type="file" id="networkfileselected" name="networkfileselected[]" multiple="multiple"/></button>
		<br>
    %else:
        
        <p><b>Samples:</b></p>
        <div id="datasets" class="multiselect" style="display: none;">
        </div>
    %end
<!--         <br>
        <button id="circvisClear">Clear circvis</button> -->

        <!-- <b>Platforms</b><br>
        {{platform if 'platform' in locals() else '-'}}<br> -->
		<br><b>Export</b><br>
		<button id="export_table">Table</button>
		<button id="export_svg">SVG</button><br>
		<br>

        
    </div>
</div>

<script type="text/javascript">

    var parsed_data = { complete_network: [], features:[], located_features:[], network:[], unlocated:[], unlocated_features:[]};


    function initMultiSelect () {
        var data = '{{data}}';
        if (data !=='') {
            $('#datasets').css('display','block')
            var data_root = '{{data}}';
            $.get(data_root + '/filelist.txt', function(data) {
                var files = data.split('\n');

                for (var i = 0; i < files.length-1; i++) {
                    if (files[i] !== 'filelist.txt'){
                        var select_str = '<a class="dataset" style="cursor: pointer;">' + files[i].replace('.tsv','') + '</a> <br>';
                        $('#datasets').append(select_str);
                    };
                };
                $('.multiselect').multiselect();
                //Set default dataset
                var default_data = '{{default}}';
                if (default_data === '') {
                    default_data = 'BPH_337';
                }
                updateCircvis($('.dataset').filter(function(){
                    return $(this).text() === default_data; 
                }));
                $('.dataset').filter(function(){
                    return $(this).text() === default_data; 
                }).addClass("multiselect-on");
            }); 
        }
    };



    /*
    File selection and secured reading are supported by only modern browsers and
    do not work on most IE versions
    */
    function UploadHandler(e) {

        var browser = navigator.appName;
        if (browser.indexOf("Internet Explorer") != -1){
            var versiontk = navigator.appVersion.split(";")[1].split(" ");
            var version = parseInt(versiontk[versiontk.length-1]);
            if (version < 10){
                console.log("Error", "Internet Explorer " + version + " does not support secured reading of local files, please use HTML5 compatible browsers such as Firefox and Chrome");
                return;
            }
            if (version >= 10) {
                console.log("Warn", "Internet Explorer " + version + " supposed to support secured reading of local files, but we have only done limited testing because of OS limitations. We highly recommend using Firefox or Chrome or Safari");
            }
        }
        var files = e.target.files || e.dataTransfer.files;
        for (var i = 0, f; f = files[i]; i++) {
            var reader = new FileReader();
            reader.onload=function(e){
                var data = this.result;
                drawCircvis(data);
            };
            reader.readAsText(f);
        }
    };

    function updateCircvis (that) {
        
        var data_root = '{{data}}';
        var datastring = "";

        var filename = that.text();
        $.get(data_root + '/'+filename+'.tsv', function(data) {
            datastring = datastring + data;
            drawCircvis(data);
        });
    };

    jQuery.fn.multiselect = function() {
        $(this).each(function() {
            var checkboxes = $(this).find('a');
            checkboxes.each(function() {
                var checkbox = $(this);
     
                // Highlight links that the user selects
                checkbox.click(function() {
                    $('.dataset').removeClass("multiselect-on");
                    checkbox.addClass("multiselect-on");
                    updateCircvis($(this));
                });
            });
        });
    };

    function clearCircvis () {
        document.getElementById("networkfileselected").value = ""
        drawCircvis("");
    };
	
    $(function() {

        // $("#circvisClear").button().click(function(e) {clearCircvis();});

        var data = '{{data}}';
        if (data===''){
        var networkfileselected = document.getElementById("networkfileselected");
        networkfileselected.addEventListener("change", UploadHandler, false); 
        }
        var data_root = '{{data}}';
        
        var data ="";
        drawCircvis(data);

        initMultiSelect();
		
		$('#export_table').button().click(function() {});
		
		$('#export_svg').button().click(function() {
			export_svg('circvis.svg') });
    });

</script>
