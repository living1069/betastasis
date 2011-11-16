Ext.ns('org.betastasis.visualizations.mutation_map');

console.log('org.betastasis.visualizations.mutation_map! available at:' + org.betastasis.visualizations.mutation_map);

var MutationMapContainer = Ext.extend(Object, {
    currentState: {},

    constructor: function(container) {
        console.log("org.betastasis.visualizations.mutation_map.MutationMapContainer(" + container + ")");
        Ext.apply(this, {contentEl:container});
    },

    logo: {
       url: "https://betastasis.googlecode.com/files/logo.png",
       label: "Betastasis Mutation map"
    },

    draw: function(data, options) {
        console.log("org.betastasis.visualizations.mutation_map.MutationMapContainer(" + data + "," + options + ")");

        var container = Ext.getDom(this.contentEl);

        Ext.Ajax.request({
            url: data.uri,
            method: "get",
            success: function(o) {
                var data = Ext.util.JSON.decode(o.responseText);
                var vis = new ProbesetVis(container);
                vis.draw(data, options);
            }
        });
    },

    GetState: function() {
        console.log("org.betastasis.visualizations.mutation_map.MutationMapContainer.GetState(): TODO: gather state from contained visualization");
        return this.currentState;
    },

    SetState: function(state) {
        this.currentState = state;
        console.log("org.betastasis.visualizations.mutation_map.MutationMapContainer.SetState(" + state + "): TODO: pass state on to contained visualization");
    }
});

org.betastasis.visualizations.mutation_map.Container = MutationMapContainer;

org.systemsbiology.pages.apis.events.MessageBus.fireEvent("publish", {
    key: "org.betastasis.visualizations.mutation_map.Container",
    payload: org.betastasis.visualizations.mutation_map.Container
});
