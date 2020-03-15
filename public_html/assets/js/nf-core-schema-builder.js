//
// nf-core-schema-builder.js
// Custom javascript for the nf-core JSON Schema Builder interface
//

// TODO - make Enter and Tab / Shift+Tab move around fields. Shift + up/down to move up and down.

// TODO - JSON Schema required array should be in parent object only? needs thought for groups

// Global variables
var schema = '';
var new_param_idx = 1;
var new_group_idx = 1;

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
    $('#schema-builder').html( generate_obj(schema['properties']['params']['properties'], 1) );
    init_group_sortable();

    // Add parameter button
    $('.add-param-btn').click(function(e){
        var new_id = 'new_param_'+new_param_idx;
        while (Object.keys(schema['properties']['params']['properties']).indexOf(new_id) != -1) {
            new_param_idx += 1;
            new_id = 'new_param_'+new_param_idx;
        }
        var new_param = {
            "type": "string",
            "description": "",
            "default": ""
        };
        schema['properties']['params']['properties'][new_id] = new_param;
        param_row = $( generate_param_row(new_id, new_param) );
        param_row.appendTo('#schema-builder').find('.param_id').select();
        scroll_to( param_row );
        new_param_idx += 1;

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    // Add group button
    $('.add-group-btn').click(function(e){
        var new_id = 'New Group '+new_group_idx;
        while (Object.keys(schema['properties']['params']['properties']).indexOf(new_id) != -1) {
            new_group_idx += 1;
            new_id = 'New Group '+new_group_idx;
        }
        var new_param = {
            "type": "object",
            "description": "",
            "default": "",
            "properties": {}
        };
        schema['properties']['params']['properties'][new_id] = new_param;
        param_row = $( generate_group_row(new_id, new_param) );
        param_row.appendTo('#schema-builder').find('.param_id').select();
        scroll_to( param_row );
        init_group_sortable();
        new_group_idx += 1;

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
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
            $('.add-param-btn').attr('disabled', true);
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
        } else {
            $('.add-param-btn').attr('disabled', false);
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
                if(!validate_id(new_id, id)){
                    $(this).val(id);
                    $(this).focus();
                } else {

                    // Do it in a slightly odd way to preserve key order
                    var new_schema = JSON.parse(JSON.stringify(schema));
                    new_schema['properties']['params']['properties'] = {};
                    for(k in schema['properties']['params']['properties']){
                        var new_k = k;
                        if(k == id){ new_k =  new_id};
                        new_schema['properties']['params']['properties'][new_k] = schema['properties']['params']['properties'][k];
                    }
                    schema = new_schema;

                    id = new_id;
                    row.data('id', id);
                    if(row.hasClass('schema_group_row')){
                        var group = row.closest('.schema_group');
                        group.data('id', id);
                        group.find('.card-body').data('id', id);
                    }
                }
            }
        }

        // Update param keys if changed
        if($(this).hasClass('param_key')){
            var param_key = $(this).data('param_key');
            var param = find_param_in_schema(id);
            var new_param = JSON.parse(JSON.stringify(param));
            new_param[param_key] = $(this).val().trim();

            // Validate
            if(!validate_param(new_param)){
                // Replace value with old copy from schema
                $(this).val( param[param_key] );
                $(this).focus();
            } else {

                param = new_param;

                // Type has changed - rebuild row
                if(param_key == 'type'){

                    // If now a boolean and defult is not string 'True', set to False
                    if($(this).val() == 'boolean'){
                        var param_default = row.find("input[data-param_key='default']").val().trim();
                        if(param_default.toLowerCase() != 'true'){
                            param['default'] = 'False';
                        }
                    }

                    // Remove special settings if not supported by the type
                    if($(this).val() != 'string'){
                        delete param['pattern'];
                    }
                    if($(this).val() != 'range'){
                        delete param['minimum'];
                        delete param['maximum'];
                    }

                    // Validate and empty default before we build the HTML to avoid warnings
                    var focus_default = false;
                    if(!validate_param(param)){
                        param['default'] = '';
                        focus_default = true;
                    }

                    row.replaceWith(generate_param_row(id, param));

                    if(focus_default){
                        $(".schema_row[data-id='"+id+"'] .param_key[data-param_key='default']").removeAttr('value').focus();
                    }
                }

                update_param_in_schema(id, param);
            }

        }

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Sorting - element has been moved
    //
    $('#schema-builder, .schema_group .card-body').on('sortstop', function(e, ui){
        // Don't actually need to know where it landed - just rebuild schema from the DOM
        var new_schema = JSON.parse(JSON.stringify(schema));
        new_schema['properties']['params']['properties'] = {};
        $('.schema_row').each(function(idx, row){
            var id = $(row).data('id');
            var param = JSON.parse(JSON.stringify(find_param_in_schema(id)));

            // Check if we are inside a group
            if ($(this).parent('.card-body').length) {
                var group_id = $(this).parent().data('id');
                new_schema['properties']['params']['properties'][group_id]['properties'][id] = param;
            } else {
                new_schema['properties']['params']['properties'][id] = param;
                // If a group, delete contents of that as well
                if(new_schema['properties']['params']['properties'][id].hasOwnProperty('properties')){
                    new_schema['properties']['params']['properties'][id]['properties'] = {};
                }
            }
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
        if(!schema['properties']['params'].hasOwnProperty('required')){
            schema['properties']['params']['required'] = [];
        }
        if(is_required){
            schema['properties']['params']['required'].push(id);
        } else {
            var idx = schema['properties']['params']['required'].indexOf(id);
            if (idx !== -1) {
                schema['properties']['params']['required'].splice(idx, 1);
            }
        }
        // Remove required array if empty
        if(schema['properties']['params']['required'].length == 0){
            delete schema['properties']['params']['required'];
        }
        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Hidden - hidden checkbox pressed
    //
    $('#schema-builder').on('change', 'input.param_hidden', function(){
        var id = $(this).closest('.schema_row').data('id');
        var is_required = $(this).is(':checked');

        // Find and update param
        var param = find_param_in_schema(id);
        if(is_required){
            param['hidden'] = true;
        } else {
            delete param['hidden'];
        }
        update_param_in_schema(id, param);

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Settings modal
    //
    $('#schema-builder').on('click', '.schema_row_config', function(){
        // Get row
        var row = $(this).closest('.schema_row');
        var id = row.data('id');
        var param = find_param_in_schema(id);

        // Build modal
        var modal_header = 'params.<span>'+id+'</span>';
        if(param['type'] == 'object'){
            modal_header = '<span>'+id+'</span>';
        }
        $('#settings_modal .modal-title').html(modal_header);
        $('#settings_enum, #settings_pattern, #settings_minimum, #settings_maximum, #settings_fa_icon').val('');
        $('.settings_enum_group, .settings_pattern_group, .settings_minmax_group').hide();

        if(['boolean', 'object'].indexOf(param['type']) == -1){
            $('.settings_enum_group').show();
        }
        if(param['type'] == 'string'){
            $('.settings_pattern_group').show();
        }
        if(param['type'] == 'range'){
            $('.settings_minmax_group').show();
        }

        // Fill modal boxes
        if(param['enum'] instanceof Array){
            $('#settings_enum').val( param['enum'].join('|') );
        }
        if(param.hasOwnProperty('pattern')){
            $('#settings_pattern').val( param['pattern'] );
        }
        if(param.hasOwnProperty('minimum')){
            $('#settings_minimum').val( param['minimum'] );
        }
        if(param.hasOwnProperty('maximum')){
            $('#settings_maximum').val( param['maximum'] );
        }
        if(param.hasOwnProperty('fa_icon')){
            $('#settings_fa_icon').val( param['fa_icon'] );
        }

        $('#settings_modal').modal('show');
    });

    //
    // Settings Modal - save button
    //
    $('#settings_save').click(function(e){

        var id = $('#settings_modal .modal-title span').text();
        var param = find_param_in_schema(id);

        var settings = {};
        settings.enum = $('#settings_enum').val().trim().split('|');
        // Trim whitespace from each element and remove empties
        $.map(settings.enum, $.trim);
        settings.enum = settings.enum.filter(function (el) { return el.length > 0; });
        settings.pattern = $('#settings_pattern').val().trim();
        settings.minimum = $('#settings_minimum').val().trim();
        settings.maximum = $('#settings_maximum').val().trim();
        settings.fa_icon = $('#settings_fa_icon').val().trim();

        // Validate inputs
        if(settings.minimum.length > 0){
            if(isNaN(parseFloat(settings.minimum))){
                alert('Error: Minimum value must be numeric');
                e.preventDefault();
                e.stopPropagation();
            }
        }
        if(settings.maximum.length > 0){
            if(isNaN(parseFloat(settings.maximum))){
                alert('Error: Maximum value must be numeric');
                e.preventDefault();
                e.stopPropagation();
            }
        }
        if(settings.minimum.length > 0 && settings.maximum.length > 0){
            if(settings.maximum < settings.minimum){
                alert('Error: Maximum value must be more than minimum');
                e.preventDefault();
                e.stopPropagation();
            }
        }
        // Check that the font-awesome icon looks right
        if(settings.fa_icon.length > 0){
            var re = new RegExp('<i class="fa[a-z -]+"><\/i>');
            if(!re.test(settings.fa_icon)){
                alert('Error: FontAwesome icon must match the regex: <i class="fa[a-z -]+"><\/i>');
                return false;
            } else {
               // Update the icon in the row
               $(".schema_row[data-id='"+id+"'] .param_fa_icon i").replaceWith($(settings.fa_icon));
           }
       } else {
           // Update the icon in the row
           $(".schema_row[data-id='"+id+"'] .param_fa_icon i").replaceWith($('<i class="far fa-question-circle fa-fw param_fa_icon_missing"></i>'));
       }

        // Update the schema
        for (var key in settings) {
            if(settings[key].length > 0){
                param[key] = settings[key];
            } else {
                delete param[key];
            }
        }

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));


    });
    // Revalidate default value once modal settings changed
    $('#settings_modal').on('hidden.bs.modal', function (e) {
        var id = $('#settings_modal .modal-title span').text();
        var param = find_param_in_schema(id);
        // It may have been deleted
        if(param){
            if(!validate_param(param)){
                $(".schema_row[data-id='"+id+"'] .param_key[data-param_key='default']").focus();
            }
        }
    });

    //
    // Settings Modal - delete button
    //
    $('#settings_delete').click(function(e){
        var id = $('#settings_modal .modal-title span').text();
        var row_el = $('.schema_row[data-id="'+id+'"]');
        var group_el = $('.schema_group[data-id="'+id+'"]');

        for(k in schema['properties']['params']['properties']){
            // Check if group
            if(schema['properties']['params']['properties'][k].hasOwnProperty('properties')){
                for(j in schema['properties']['params']['properties'][k]['properties']){
                    if(j == id){
                        delete schema['properties']['params']['properties'][k]['properties'][j];
                    }
                    // Move contents of the group out if we're going to delete it
                    if(k == id){
                        schema['properties']['params']['properties'][j] = schema['properties']['params']['properties'][k]['properties'][j];
                        // NOT WORKING
                        $('.schema_row[data-id="'+j+'"]').insertBefore(group_el);
                    }
                }
            }
            if(k == id){
                delete schema['properties']['params']['properties'][k];
            }
        }

        row_el.remove();
        group_el.remove();

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Copy schema button
    //
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
    $([document.documentElement, document.body]).animate({
        scrollTop: el_offset
    }, 500);
}

function generate_obj(obj, level){
    var results = '';
    for (var id in obj){
        if (obj.hasOwnProperty(id)) {

            // Groups
            if(obj[id]['type'] == 'object'){

                // Generate child rows
                var child_params = '';
                for (var child_id in obj[id]['properties']){
                    if (obj[id]['properties'].hasOwnProperty(child_id)) {
                        child_params += generate_param_row(child_id, obj[id]['properties'][child_id]);
                    }
                }
                results += generate_group_row(id, obj[id], child_params);
            }

            // Regular rows
            else {
                results += generate_param_row(id, obj[id]);
            }
        }
    }
    return results;
}

function generate_param_row(id, param){

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
        if(/^-?[\d\.]+$/.test(param['minimum'])){
            attrs += ' min="'+param['minimum']+'"';
        }
        if(/^-?[\d\.]+$/.test(param['maximum'])){
            attrs += ' max="'+param['maximum']+'"';
        }
        if(param['default'] != undefined){
            attrs += ' value="'+param['default']+'"';
        }
        default_input = '<input '+attrs+' class="param_key" data-param_key="default">';
    }

    var is_required = false;
    // Check that the required array exists
    if(schema['properties']['params'].hasOwnProperty('required')){
        if (schema['properties']['params']['required'].indexOf(id) !== -1) {
            is_required = true;
        }
    }

    var is_hidden = false;
    if(param['hidden']){
        is_hidden = true;
    }

    var fa_icon = '<i class="far fa-question-circle fa-fw param_fa_icon_missing"></i>';
    if(param['fa_icon'] != undefined){
        fa_icon = $(param['fa_icon']).addClass('fa-fw').get(0).outerHTML;
    }


    var results = `
    <div class="row schema_row" data-id="`+id+`">
        <div class="col-sm-auto align-self-center d-none d-sm-block schema_row_grabber">
            <i class="fas fa-grip-vertical"></i>
        </div>
        <div class="col-sm-auto align-self-center d-none d-sm-block param_fa_icon ">`+fa_icon+`</div>
        <div class="col schema-id">
            <label>ID
                <input type="text" class="text-monospace param_id" value="`+id+`">
            </label>
        </div>
        <div class="col-sm-3">
            <label>Description
                <input type="text" class="param_key" data-param_key="description" value="`+param['description']+`">
            </label>
        </div>
        <div class="col-sm-1">
            <label>Type
                <select class="param_key" data-param_key="type">
                    <option `+(param['type'] == 'string' ? 'selected="selected"' : '')+` value="string">string</option>
                    <option `+(param['type'] == 'number' ? 'selected="selected"' : '')+` value="number">number</option>
                    <option `+(param['type'] == 'integer' ? 'selected="selected"' : '')+` value="integer">integer</option>
                    <option `+(param['type'] == 'range' ? 'selected="selected"' : '')+` value="range">range</option>
                    <option `+(param['type'] == 'boolean' ? 'selected="selected"' : '')+` value="boolean">boolean</option>
                </select>
            </label>
        </div>
        <div class="col-sm-3">
            <label>Default `+default_input+`</label>
        </div>
        <div class="col-sm-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>Required
                <input type="checkbox" `+(is_required ? 'checked="checked"' : '')+`" class="param_required">
            </label>
            `)+`
        </div>
        <div class="col-sm-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>Hidden
                <input type="checkbox" `+(is_hidden ? 'checked="checked"' : '')+`" class="param_hidden">
            </label>
            `)+`
        </div>
        <div class="col-sm-auto align-self-center schema_row_config">
            <i class="fas fa-cog"></i>
        </div>
    </div>`;

    return results;
}

function generate_group_row(id, param, child_params){

    if(child_params == undefined){
        child_params = '';
    }

    var results = `
    <div class="card schema_group" data-id="`+id+`">
        <div class="card-header p-0">
            <div class="row schema_row schema_group_row mb-0" data-id="`+id+`">
                <div class="col-sm-auto align-self-center schema_row_grabber d-none d-sm-block">
                    <i class="fas fa-grip-vertical"></i>
                </div>
                <div class="col schema-id">
                    <label>Title
                        <input type="text" class="text-monospace param_id" value="`+id+`">
                    </label>
                </div>
                <div class="col">
                    <label>Description
                        <input type="text" class="param_key" data-param_key="description" value="`+param['description']+`">
                    </label>
                </div>
                <div class="col-sm-auto align-self-center schema_row_config">
                    <i class="fas fa-cog"></i>
                </div>
            </div>
        </div>
        <div class="card-body" data-id="`+id+`">`+child_params+`</div>
    </div>
    `;

    return results;
}

function init_group_sortable(){
    // Initialise sortable functionality
    // Done in a function as these can be dynamically created, so may need to be initiliased more than once

    // Main body
    $('#schema-builder').sortable({
        handle: '.schema_row_grabber',
        tolerance: 'pointer',
        placeholder: 'schema_row_move_placeholder alert alert-warning',
        connectWith: '.schema_group .card-body'
    });

    // Object Groups
    $(".schema_group .card-body").sortable({
        handle: '.schema_row_grabber',
        tolerance: 'pointer',
        placeholder: 'schema_row_move_placeholder alert alert-warning',
        connectWith: '#schema-builder, .schema_group .card-body'
    });

    // Listeners to prevent nested groups
    $(".schema_group .card-body").on("sortreceive", function(e, ui) {
        if(ui.item.hasClass('schema_group')){
            ui.sender.sortable('cancel');
        }
    });
}

function validate_id(id, old_id){

    var param = false;
    var is_object = false;

    // Get param if we have the old ID
    if(old_id !== undefined){
        param = find_param_in_schema(old_id);
        is_object = (param['type'] == 'object');
    }

    // Check that the ID is simple
    if(!is_object){
        var re = new RegExp("^[a-zA-Z0-9_]+$");
        if(!re.test(id)){
            alert('Error: Parameter ID must be just alphanumeric / underscores');
            return false;
        }
    }

    // Check that the ID is not a duplicate
    var num_hits = 0;
    // Simple case - not in a group
    if(schema['properties']['params']['properties'].hasOwnProperty(id)){
        num_hits += 1;
    }
    // Iterate through groups, looking for ID
    for(k in schema['properties']['params']['properties']){
        if(schema['properties']['params']['properties'][k].hasOwnProperty('properties')){
            if(schema['properties']['params']['properties'][k]['properties'].hasOwnProperty(id)){
                num_hits += 1;
            }
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
        if(param.hasOwnProperty('minimum') && !isNaN(parseFloat(param['minimum']))){
            if(parseFloat(param['default']) < parseFloat(param['minimum'])){
                alert('Error: Default value "'+param['default']+'" must be greater than or equal to minimum value: '+param['minimum']);
                return false;
            }
        }
        if(param.hasOwnProperty('maximum') && !isNaN(parseFloat(param['maximum']))){
            if(parseFloat(param['default']) > parseFloat(param['maximum'])){
                alert('Error: Default value "'+param['default']+'" must be less than or equal to maximum value: '+param['maximum']);
                return false;
            }
        }
    }

    // Empty defaults are always ok
    if(param.hasOwnProperty('default') && param['default'].length == 0){
        return true;
    }

    // Check that numbers and ranges are numbers
    if(['number', 'range'].includes(param['type'])){
        var default_float = parseFloat(param['default']);
        if(String(default_float) !== String(param['default'])){
            alert('Error: Default value "'+param['default']+'" is not a number');
            return false;
        }
    }

    // Check that integers are integers
    if(param['type'] == 'integer'){
        var default_int = parseInt(param['default']);
        if(String(default_int) !== String(param['default'])){
            alert('Error: Default value "'+param['default']+'" is not an integer');
            return false;
        }
    }

    // Check that default matches enum
    if(param['enum'] instanceof Array){
        if(param['enum'].indexOf(param['default']) == -1){
            alert('Error: Default value "'+param['default']+'" must be one of the Enumerated values: '+param['enum'].join(', '));
            return false;
        }
    }

    // Check that default matches regex pattern
    if(param.hasOwnProperty('pattern')){
        var re = new RegExp(param['pattern']);
        if(!re.test(param['default'])){
            alert('Error: Default value "'+param['default']+'" must match the parameter pattern regex: '+param['pattern']);
            return false;
        }
    }

    return true;
}


function find_param_in_schema(id){
    // Given an ID, find the param schema even if it's in a group
    // Assumes max one level of nesting and unique IDs everywhere

    // Simple case - not in a group
    if(schema['properties']['params']['properties'].hasOwnProperty(id)){
        return schema['properties']['params']['properties'][id];
    }

    // Iterate through groups, looking for ID
    for(k in schema['properties']['params']['properties']){
        // Check if group
        if(schema['properties']['params']['properties'][k].hasOwnProperty('properties')){
            if(schema['properties']['params']['properties'][k]['properties'].hasOwnProperty(id)){
                return schema['properties']['params']['properties'][k]['properties'][id];
            }
        }
    }

    console.warn("Could not find param '"+id+"'");
}

function update_param_in_schema(id, new_param){
    // Given an ID, find the param schema even if it's in a group
    // Assumes max one level of nesting and unique IDs everywhere

    // Simple case - not in a group
    if(schema['properties']['params']['properties'].hasOwnProperty(id)){
        schema['properties']['params']['properties'][id] = new_param;
        return true;
    }

    // Iterate through groups, looking for ID
    for(k in schema['properties']['params']['properties']){
        // Check if group
        if(schema['properties']['params']['properties'][k].hasOwnProperty('properties')){
            if(schema['properties']['params']['properties'][k]['properties'].hasOwnProperty(id)){
                schema['properties']['params']['properties'][k]['properties'][id] = new_param;
                return true;
            }
        }
    }

    console.warn("Could not find param '"+id+"'");
}
