//
// nf-core-schema-builder.js
// Custom javascript for the nf-core JSON Schema Builder interface
//

// Global variables
var schema = '';

$(function () {

    // Show the cache expiration time in local timezone
    $('.cache_expires_at span').text(
        moment.unix($('.cache_expires_at span').text()).calendar().replace(/^\w/, function (chr) { return chr.toLowerCase(); })
    );
    $('.cache_expires_at').show();

    // Parse initial JSON Schema
    schema = JSON.parse($('#json_schema').text());

    // Build the schema builder
    $('.schema-builder').html( generate_obj(schema['properties']['input']['properties'], 1) );

    // Listeners to update on change
    $('.schema-builder').on('change', 'input, select', function(){
        var row = $(this).closest('.schema_row');

        // Parse data attributes
        var id = row.data('id');

        // Update ID if changed
        if($(this).hasClass('param_id')){
            var new_id = $(this).val();

            // TODO - doesn't handle objects / groups
            // Do it in a slightly odd way to preserver key order
            var new_schema = JSON.parse(JSON.stringify(schema));
            new_schema['properties']['input']['properties'] = {};
            for(k in schema['properties']['input']['properties']){
                var new_k = k;
                if(k == id){ new_k =  new_id};
                new_schema['properties']['input']['properties'][new_k] = schema['properties']['input']['properties'][k];
            }
            schema = new_schema;

            id = new_id;
            row.data('id', id);
        }

        // Update param keys if changed
        if($(this).hasClass('param_key')){
            var param_key = $(this).data('param_key');
            schema['properties']['input']['properties'][id][param_key] = $(this).val();

            // Type has changed - rebuild row
            if(param_key == 'type'){

                // If now a boolean and defult is not string 'True', set to False
                if($(this).val() == 'boolean'){
                    var param_default = row.find("input[data-param_key='default']").val();
                    if(param_default.toLowerCase() != 'true'){
                        schema['properties']['input']['properties'][id]['default'] = 'False';
                    }
                }

                row.replaceWith(generate_row(id, schema['properties']['input']['properties'][id]));
            }

        }

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

});


function generate_obj(obj, level){
    var results = '';
    for (var id in obj){
        if (obj.hasOwnProperty(id)) {
            results += generate_row(id, obj[id]);
        }
    }
    return results;
}

function generate_row(id, param){

    results = `
    <div class="row schema_row" data-id="`+id+`">
        <div class="col-sm-auto align-self-center schema_row_grabber">
            <i class="fas fa-grip-vertical"></i>
        </div>
        <div class="col-sm-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>Required
                <input type="checkbox" checked="`+param['default']+`">
            </label>
            `)+`
        </div>
        <div class="col schema-id">
            <label>ID
                <input type="text" class="text-monospace param_id" value="`+id+`">
            </label>
        </div>
        <div class="col-sm-2">
            <label>Type
                <select class="param_key" data-param_key="type">
                    <option `+(param['type'] == 'string' ? 'selected="selected"' : '')+` value="string">string</option>
                    <option `+(param['type'] == 'number' ? 'selected="selected"' : '')+` value="number">number</option>
                    <option `+(param['type'] == 'boolean' ? 'selected="selected"' : '')+` value="boolean">boolean</option>
                    <option `+(param['type'] == 'object' ? 'selected="selected"' : '')+` value="object">object (group)</option>
                </select>
            </label>
        </div>
        <div class="col-sm-3">
            <label>Description
                <input type="text" class="param_key" data-param_key="description" value="`+param['description']+`">
            </label>
        </div>
        <div class="col-sm-3">
            <label>Default
                `+(param['type'] == 'boolean' ? `
                <select class="param_key" data-param_key="default">
                    <option `+(param['default'] == 'True' ? 'selected="selected"' : '')+`>True</option>
                    <option `+(param['default'] == 'False' ? 'selected="selected"' : '')+`>False</option>
                </select>
                ` : '<input type="text" class="param_key" data-param_key="default" value="'+param['default']+'">')+`
            </label>
        </div>
        <div class="col-sm-auto align-self-center schema_row_config">
            <i class="fas fa-cog"></i>
        </div>
    </div>`;

    return results;
}

function validate_param(){

}
