function object(d){var s=function(){};s.prototype=d;return new s();}var QuickSelect;(function($){QuickSelect=function(d,f){var self=this;d=$(d);d.attr('autocomplete','off');self.options=f;self.AllItems={};var g=false,h=-1,j=false,k,l,m,n=false,o,p;if(/MSIE (\d+\.\d+);/.test(navigator.userAgent)){if(Number(RegExp.$1)<=7)n=true;}o=$('<div class="'+f.resultsClass+'" style="display:block;position:absolute;z-index:9999;"></div>').hide();p=$('<iframe />');p.css({border:'none',position:'absolute'});if(f.width>0){o.css("width",f.width);p.css("width",f.width);}$('body').append(o);o.hide();if(n)$('body').append(p);self.getLabel=function(A){return A.label||(typeof(A)==='string'?A:A[0])||'';};var r=function(A){return A.values||(A.value?[A.value]:(typeof(A)==='string'?[A]:A))||[];};var t=function(A){var B=$('li',o);if(!B)return;if(typeof(A)==="number")h=h+A;else h=B.index(A);if(h<0)h=0;else if(h>=B.size())h=B.size()-1;B.removeClass(f.selectedClass);$(B[h]).addClass(f.selectedClass);if(f.autoFill&&self.last_keyCode!=8){d.val(l+$(B[h]).text().substring(l.length));var C=l.length,D=d.val().length,E=d.get(0);if(E.createTextRange){var F=E.createTextRange();F.collapse(true);F.moveStart("character",C);F.moveEnd("character",D);F.select();}else if(E.setSelectionRange){E.setSelectionRange(C,D);}else if(E.selectionStart){E.selectionStart=C;E.selectionEnd=D;}E.focus();}};var u=function(){if(m){clearTimeout(m);}d.removeClass(f.loadingClass);if(o.is(":visible"))o.hide();if(p.is(":visible"))p.hide();h=-1;};self.selectItem=function(A,B){if(!A){A=document.createElement("li");A.item='';}var C=self.getLabel(A.item),D=r(A.item);d.lastSelected=C;d.val(C);l=C;o.empty();$(f.additionalFields).each(function(i,E){$(E).val(D[i+1]);});if(!B)u();if(f.onItemSelect)setTimeout(function(){f.onItemSelect(A);},1);return true;};var v=function(){var A=$("li."+f.selectedClass,o).get(0);if(A){return self.selectItem(A);}else{if(f.exactMatch){d.val('');$(f.additionalFields).each(function(i,B){$(B).val('');});}return false;}};var w=function(A){o.empty();if(!j||A===null||A.length===0)return u();var B=document.createElement("ul"),C=A.length,D=function(){t(this);},E=function(){},F=function(e){e.preventDefault();e.stopPropagation();self.selectItem(this);};o.append(B);if(f.maxVisibleItems>0&&f.maxVisibleItems<C)C=f.maxVisibleItems;for(var i=0;i<C;i++){var G=A[i],H=document.createElement("li");o.append(H);$(H).text(f.formatItem?f.formatItem(G,i,C):self.getLabel(G));H.item=G;if(G.className)H.className=G.className;B.appendChild(H);$(H).hover(D,E).click(F);}d.removeClass(f.loadingClass);return true;};var x=function(q,A){f.finderFunction.apply(self,[q,function(B){w(f.matchMethod.apply(self,[q,B]));A();}]);};var y=function(){var A=d.offset(),B=(f.width>0?f.width:d.width()),C=$('li',o);o.css({width:parseInt(B,10)+"px",top:A.top+d.height()+5+"px",left:A.left+"px"});if(n){p.css({width:parseInt(B,10)-2+"px",top:A.top+d.height()+6+"px",left:A.left+1+"px",height:o.height()-2+'px'}).show();}o.show();if(f.autoSelectFirst||(f.selectSingleMatch&&C.length==1))t(C.get(0));};var z=function(){if(k>=9&&k<=45){return;}var q=d.val();if(q==l)return;l=q;if(q.length>=f.minChars){d.addClass(f.loadingClass);x(q,y);}else{if(q.length===0&&(f.onBlank?f.onBlank():true))$(f.additionalFields).each(function(i,A){A.value='';});d.removeClass(f.loadingClass);o.hide();p.hide();}};o.mousedown(function(e){if(e.srcElement)g=e.srcElement.tagName!='DIV';});d.keydown(function(e){k=e.keyCode;switch(e.keyCode){case 38:e.preventDefault();t(-1);break;case 40:e.preventDefault();if(!o.is(":visible")){y();t(0);}else{t(1);}break;case 13:if(v()){e.preventDefault();d.select();}break;case 9:break;case 27:if(h>-1&&f.exactMatch&&d.val()!=$($('li',o).get(h)).text()){h=-1;}$('li',o).removeClass(f.selectedClass);u();e.preventDefault();break;default:if(m){clearTimeout(m);}m=setTimeout(z,f.delay);break;}}).focus(function(){j=true;}).blur(function(e){if(h>-1){v();}j=false;if(m){clearTimeout(m);}m=setTimeout(function(){u();if(f.exactMatch&&d.val()!=d.lastSelected){self.selectItem(null,true);}},200);});};QuickSelect.matchers={quicksilver:function(q,d){var f,g,self=this;f=(self.options.matchCase?q:q.toLowerCase());self.AllItems[f]=[];for(var i=0;i<d.length;i++){g=(self.options.matchCase?self.getLabel(d[i]):self.getLabel(d[i]).toLowerCase());if(g.score(f)>0){self.AllItems[f].push(d[i]);}}return self.AllItems[f].sort(function(a,b){a=(self.options.matchCase?self.getLabel(a):self.getLabel(a).toLowerCase());b=(self.options.matchCase?self.getLabel(b):self.getLabel(b).toLowerCase());a=a.score(f);b=b.score(f);return(a>b?-1:(b>a?1:0));});},contains:function(q,d){var f,g,self=this;f=(self.options.matchCase?q:q.toLowerCase());self.AllItems[f]=[];for(var i=0;i<d.length;i++){g=(self.options.matchCase?self.getLabel(d[i]):self.getLabel(d[i]).toLowerCase());if(g.indexOf(f)>-1){self.AllItems[f].push(d[i]);}}return self.AllItems[f].sort(function(a,b){a=(self.options.matchCase?self.getLabel(a):self.getLabel(a).toLowerCase());b=(self.options.matchCase?self.getLabel(b):self.getLabel(b).toLowerCase());var h=a.indexOf(f);var j=b.indexOf(f);return(h>j?-1:(h<j?1:(a>b?-1:(b>a?1:0))));});},startsWith:function(q,d){var f,g,self=this;f=(self.options.matchCase?q:q.toLowerCase());self.AllItems[f]=[];for(var i=0;i<d.length;i++){g=(self.options.matchCase?self.getLabel(d[i]):self.getLabel(d[i]).toLowerCase());if(g.indexOf(f)===0){self.AllItems[f].push(d[i]);}}return self.AllItems[f].sort(function(a,b){a=(self.options.matchCase?self.getLabel(a):self.getLabel(a).toLowerCase());b=(self.options.matchCase?self.getLabel(b):self.getLabel(b).toLowerCase());return(a>b?-1:(b>a?1:0));});}};QuickSelect.finders={data:function(q,d){d(this.options.data);},ajax:function(q,d){var f=this.options.ajax+"?q="+encodeURI(q);for(var i in this.options.ajaxParams){if(this.options.ajaxParams.hasOwnProperty(i)){f+="\x26"+i+"\x3d"+encodeURI(this.options.ajaxParams[i]);}}$.getJSON(f,d);}};$.fn.quickselect=function(d,f){if(d=='instance'&&$(this).data('quickselect'))return $(this).data('quickselect');d=d||{};d.data=(typeof(d.data)==="object"&&d.data.constructor==Array)?d.data:undefined;d.ajaxParams=d.ajaxParams||{};d.delay=d.delay||400;if(!d.delay)d.delay=(!d.ajax?400:10);d.minChars=d.minChars||1;d.cssFlavor=d.cssFlavor||'quickselect';d.inputClass=d.inputClass||d.cssFlavor+"_input";d.loadingClass=d.loadingClass||d.cssFlavor+"_loading";d.resultsClass=d.resultsClass||d.cssFlavor+"_results";d.selectedClass=d.selectedClass||d.cssFlavor+"_selected";d.finderFunction=d.finderFunction||QuickSelect.finders[!d.data?'ajax':'data'];if(d.finderFunction==='data'||d.finderFunction==='ajax')d.finderFunction=QuickSelect.finders[d.finderFunction];d.matchMethod=d.matchMethod||QuickSelect.matchers[(typeof(''.score)==='function'&&'\x6c'.score('\x6c')==1?'quicksilver':'contains')];if(d.matchMethod==='quicksilver'||d.matchMethod==='contains'||d.matchMethod==='startsWith')d.matchMethod=QuickSelect.matchers[d.matchMethod];if(d.matchCase===undefined)d.matchCase=false;if(d.exactMatch===undefined)d.exactMatch=false;if(d.autoSelectFirst===undefined)d.autoSelectFirst=true;if(d.selectSingleMatch===undefined)d.selectSingleMatch=true;if(d.additionalFields===undefined)d.additionalFields=$('nothing');d.maxVisibleItems=d.maxVisibleItems||-1;if(d.autoFill===undefined||d.matchMethod!='startsWith'){d.autoFill=false;}d.width=parseInt(d.width,10)||0;return this.each(function(){var g=this,h=object(d);if(g.tagName=='INPUT'){var j=new QuickSelect(g,h);$(g).data('quickselect',j);}else if(g.tagName=='SELECT'){h.delay=h.delay||10;h.finderFunction='data';var name=g.name,k=g.id,l=g.className,m=$(g).attr('accesskey'),n=$(g).attr('tabindex'),o=$("option:selected",g).get(0);h.data=[];$('option',g).each(function(i,t){h.data.push({label:$(t).text(),values:[t.value,t.value],className:t.className});});var p=$("<input type='text' class='"+l+"' id='"+k+"_quickselect' accesskey='"+m+"' tabindex='"+n+"' />");if(o){p.val($(o).text());}var r=$("<input type='hidden' id='"+k+"' name='"+g.name+"' />");if(o){r.val(o.value);}h.additionalFields=r;$(g).after(p).after(r).remove();p.quickselect(h);}});};})(jQuery);