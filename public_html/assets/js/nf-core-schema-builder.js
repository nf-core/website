//
// nf-core-schema-builder.js
// Custom javascript for the nf-core JSON Schema Builder interface
//

// TODO - handle objects / groups

// Global variables
var schema = '';

$(function () {

    // Show the cache expiration time in local timezone
    $('.cache_expires_at span').text(
        moment.unix($('.cache_expires_at span').text()).calendar().replace(/^\w/, function (chr) { return chr.toLowerCase(); })
    );
    $('.cache_expires_at').show();

    // Parse initial JSON Schema
    try {
        schema = JSON.parse($('#json_schema').text());
    } catch(e){
        alert("Error - Schema JSON could not be parsed. See the browser console for details.");
        console.log(e);
    }

    // Build the schema builder
    $('#schema-builder').html( generate_obj(schema['properties']['input']['properties'], 1) );

    // Make the rows sortable
    $('#schema-builder').sortable({
        handle: '.schema_row_grabber',
        placeholder: 'schema_row_move_placeholder alert alert-warning'
    });

    //
    // FINISHED button
    //
    // Toggle between panels
    $('.schema-panel-btn').click(function(){
        var target = $( $(this).data('target') );
        if(target.is(':hidden')){
            $('.schema-panel:visible').fadeOut('fast', function(){
                target.fadeIn('fast');
                scroll_to(target);
            });
        } else {
            // Already visible, just scroll to top
            scroll_to(target);
        }

        // Post the results to PHP when finished
        if($(this).data('target') == '#schema-finished'){
            $('#schema-send-status').text("Saving schema..");

            post_data = {
                'post_content': 'json_schema',
                'version': 'web_builder',
                'status': 'web_builder_edited',
                'api': 'true',
                'cache_id': $('#schema_cache_id').text(),
                'schema': JSON.stringify(schema)
            };
            $.post( "json_schema_build", post_data).done(function( returned_data ) {
                console.log("Sent schema to API. Response:", returned_data);
                if(returned_data.status == 'recieved'){
                    $('#schema-send-status').text("Ok, that's it - done!");
                } else {
                    $('#schema-send-status').text("Oops, something went wrong!");
                }
            });
        }
    });

    //
    // LISTENERS
    //
    // Listeners to update on change
    $('#schema-builder').on('change', 'input, select', function(){
        var row = $(this).closest('.schema_row');

        // Parse data attributes
        var id = row.data('id');

        // Update ID if changed
        if($(this).hasClass('param_id')){
            var new_id = $(this).val().trim();

            // Check if it actually changed
            if(new_id != id){
                if(!validate_id(new_id)){
                    $(this).val(id);
                    $(this).focus();
                } else {

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
            }
        }

        // Update param keys if changed
        if($(this).hasClass('param_key')){
            var param_key = $(this).data('param_key');
            var param = JSON.parse(JSON.stringify(schema['properties']['input']['properties'][id]));
            param[param_key] = $(this).val().trim();

            // Validate
            if(!validate_param(param)){
                console.log("Parameter didn't validate, resetting");
                $(this).val( schema['properties']['input']['properties'][id][param_key] );
                $(this).focus();
            } else {

                schema['properties']['input']['properties'][id] = param;

                // Type has changed - rebuild row
                if(param_key == 'type'){

                    // If now a boolean and defult is not string 'True', set to False
                    if($(this).val() == 'boolean'){
                        var param_default = row.find("input[data-param_key='default']").val().trim();
                        if(param_default.toLowerCase() != 'true'){
                            schema['properties']['input']['properties'][id]['default'] = 'False';
                        }
                    }

                    row.replaceWith(generate_row(id, schema['properties']['input']['properties'][id]));
                }
            }

        }

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Sorting - element has been moved
    //
    $('#schema-builder').on('sortstop', function(e, ui){
        // Don't actually need to know where it landed - just rebuild schema from the DOM
        var new_schema = JSON.parse(JSON.stringify(schema));
        new_schema['properties']['input']['properties'] = {};
        $('.schema_row').each(function(idx, row){
            var id = $(row).data('id');
            new_schema['properties']['input']['properties'][id] = schema['properties']['input']['properties'][id];
        });
        schema = new_schema;
        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Required - required checkbox pressed
    //
    $('#schema-builder').on('change', 'input.param_required', function(){
        var row = $(this).closest('.schema_row');
        var id = row.data('id');
        var is_required = $(this).is(':checked');
        // Check that the required array exists
        if(!schema['properties']['input'].hasOwnProperty('required')){
            schema['properties']['input']['required'] = [];
        }
        if(is_required){
            schema['properties']['input']['required'].push(id);
        } else {
            var idx = schema['properties']['input']['required'].indexOf(id);
            if (idx !== -1) {
                schema['properties']['input']['required'].splice(idx, 1);
            }
        }
        // Remove required array if empty
        if(schema['properties']['input']['required'].length == 0){
            delete schema['properties']['input']['required'];
        }
        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    // Copy schema button
    $('.copy-schema-btn').click(function(){
        // select the content
        var target = $('#json_schema');
        var currentFocus = document.activeElement;
        target.attr('disabled', false);
        target.focus();
        target.select();

        // copy the selection
        try {
            document.execCommand("copy");
        } catch(e) {
            alert('Copy action did not work - please copy schema manually')
            console.log(e);
        }
        // restore original focus
        if (currentFocus && typeof currentFocus.focus === "function") {
            currentFocus.focus();
        }
        target.attr('disabled', true);
    });

});

function scroll_to(target_el){
    var el_offset = target_el.offset().top - 124;
    var doc_offset = $(document).scrollTop();
    if(doc_offset > el_offset){
        $([document.documentElement, document.body]).animate({
            scrollTop: el_offset
        }, 500);
    }
}

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

    var default_input = '';
    if(param['type'] == 'boolean'){
        default_input = `
            <select class="param_key" data-param_key="default">
                <option `+(param['default'] == 'True' ? 'selected="selected"' : '')+`>True</option>
                <option `+(param['default'] == 'False' ? 'selected="selected"' : '')+`>False</option>
            </select>`;
    }
    if(['string', 'integer', 'number', 'range'].includes(param['type'])){
        var attrs = '';
        if(param['type'] == 'string'){
            attrs = 'type="text"';
        } else {
            attrs = 'type="number"';
        }
        if(/^[\d.]+$/.test(param['minimum'])){
            attrs += ' min="'+param['minimum']+'"';
        }
        if(/^[\d.]+$/.test(param['maximum'])){
            attrs += ' max="'+param['maximum']+'"';
        }
        if(param['default'] != undefined){
            attrs += ' value="'+param['default']+'"';
        }
        default_input = '<input '+attrs+' class="param_key" data-param_key="default">';
    }

    var is_required = false;
    if (schema['properties']['input']['required'].indexOf(id) !== -1) {
        is_required = true;
    }


    var results = `
    <div class="row schema_row" data-id="`+id+`">
        <div class="col-sm-auto align-self-center schema_row_grabber d-none d-sm-block">
            <i class="fas fa-grip-vertical"></i>
        </div>
        <div class="col-sm-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>Required
                <input type="checkbox" `+(is_required ? 'checked="checked"' : '')+`" class="param_required">
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
                    <option `+(param['type'] == 'integer' ? 'selected="selected"' : '')+` value="integer">integer</option>
                    <option `+(param['type'] == 'range' ? 'selected="selected"' : '')+` value="range">range</option>
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
            <label>Default `+default_input+`</label>
        </div>
        <div class="col-sm-auto align-self-center schema_row_config">
            <i class="fas fa-cog"></i>
        </div>
    </div>`;

    return results;
}

function validate_id(id){
    // Check that the ID is simple
    var re = new RegExp("^[a-zA-Z0-9_]+$");
    if(!re.test(id)){
        alert('Error: Parameter ID must be just alphanumeric / underscores');
        return false;
    }

    // Check that the ID is not a duplicate
    // TODO - ignores groups
    var num_hits = 0;
    for(k in schema['properties']['input']['properties']){
        if(k == id){
            num_hits += 1;
        }
    }
    if(num_hits > 0){
        alert('Error: New parameter ID is a duplicate');
        return false;
    }

    return true;
}

function validate_param(param){

    // Check that the minimum and maximum is valid
    if(['integer', 'number', 'range'].includes(param['type'])){
        if(/^[\d.]+$/.test(param['minimum'])){
            if(param['default'] < param['minimum']){
                alert('Error: Value must be greater than or equal to '+param['minimum']);
                return false;
            }
        }
        if(/^[\d.]+$/.test(param['maximum'])){
            if(param['default'] > param['maximum']){
                alert('Error: Value must be less than or equal to '+param['minimum']);
                return false;
            }
        }
    }

    // Check integers are integers
    if(param['type'] == 'integer'){
        var default_int = parseInt(param['default']);
        if(String(default_int) !== String(param['default'])){
            alert('Error: Value is not an integer');
            return false;
        }
    }

    return true;
}
