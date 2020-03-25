//
// nf-core-schema-builder.js
// Custom javascript for the nf-core JSON Schema Builder interface
//

// Global variables
var schema = '';
var new_param_idx = 1;
var new_group_idx = 1;
var help_text_icon_template = '<i class="fas fa-comment-slash help_text_icon help_text_icon_no_text" data-toggle="tooltip" data-html="true" data-placement="right" data-delay="500" title="Does not have any help text"></i>';
var no_help_text_icon = '<i class="fas fa-comment-dots help_text_icon" data-toggle="tooltip" data-html="true" data-placement="right" data-delay="500" title="Has help text"></i>';
var prev_focus = false;

$(function () {

    // Show the cache expiration time in local timezone
    $('.cache_expires_at span').text(
        moment.unix($('.cache_expires_at span').text()).calendar().replace(/^\w/, function (chr) { return chr.toLowerCase(); })
    );
    $('.cache_expires_at').show();

    //
    // FontAwesome icon picker
    //
    // Get the list of font-awesome icons
    $.getJSON('assets/js/fa-icons.json', function(fa_icons) {

        // Make the popover
        var popover_template = `
            <div class="popover fa_icon_picker" role="tooltip">
                <div class="arrow"></div>
                <div class="popover-header"></div>
                <div class="popover-body"></div>
            </div>
        `;
        var popover_title = '<input type="text" class="form-control fa_icon_picker_input" placeholder="Search">';
        var popover_content = '<div class="d-flex flex-wrap">';
        for(icon in fa_icons){
            popover_content += '<button class="btn btn-sm btn-light m-1" data-searchterms="'+fa_icons[icon]+'" data-classname="'+icon+'"><i class="'+icon+' fa-2x fa-fw"></i></button>';
        }
        popover_content += '</div>';

        // Initialise Bootstrap popover for icon picker
        $('#schema-builder').popover({
            selector: '.param_fa_icon',
            html: true,
            sanitize: false,
            placement: 'right',
            offset: 100,
            animation: false,
            template: popover_template,
            title: popover_title,
            content: popover_content
        });

        // Listener for when the popover is triggered
        // Needs selector class instead of root class.
        $('.param_fa_icon').on('show.bs.popover', function () {
            // Only show one popover at a time
            $('.param_fa_icon').popover('hide');
            // Reset the selected icon button classes
            $('.fa_icon_picker .popover-body .btn').removeClass('btn-success').addClass('btn-light');
        });

        // Focus the search bar when triggered
        $('.param_fa_icon').on('shown.bs.popover', function () {
            var row = $(this).closest('.schema_row');
            var id = row.data('id');
            var param = find_param_in_schema(id);
            prev_focus = $(this);
            // Highlight the picked icon
            if(param.hasOwnProperty('fa_icon') && param.fa_icon.length > 0){
                var matching_btn = $('.fa_icon_picker .popover-body .btn[data-classname="'+param.fa_icon+'"]').removeClass('btn-light').addClass('btn-success');
            }
            $('.fa_icon_picker:visible').data('id', id);
            // Focus the input once we've finished showing the popover
            setTimeout(function(){
                $('.fa_icon_picker_input').focus();
            }, 10);
        });

        // Filter icons
        $('body').on('keyup', '.fa_icon_picker_input', function(e){
            var search_term = $(this).val().trim();
            var popover = $(this).closest('.fa_icon_picker');
            if(search_term.length == 0){
                popover.find(".popover-body div button").show();
            } else {
                popover.find(".popover-body div button").hide();
                popover.find(".popover-body div button[data-searchterms*='"+search_term+"'], .popover-body div button[data-classname*='"+search_term+"']").show();
            }
        });

        // Icon clicked
        $('body').on('click', '.fa_icon_picker button', function(e){
            var class_name = $(this).data('classname');

            // Update schema
            var id = $('.fa_icon_picker:visible').data('id');
            var param = find_param_in_schema(id);
            param.fa_icon = class_name;
            update_param_in_schema(id, param);

            // Update printed schema in page
            $('#json_schema').text(JSON.stringify(schema, null, 4));

            // Update form
            $('.schema_row[data-id="'+id+'"] .param_fa_icon i').removeClass().addClass(class_name+' fa-fw');
            $('.param_fa_icon').popover('hide');
            prev_focus.focus();
        });

        // Dismiss popover if click elsewhere
        $(document).click(function(e) {
            var target = $(e.target);
            // Is a click outside the icon popover and button
            if(!target.closest('.fa_icon_picker').length && !target.closest('.param_fa_icon').length){
                // Do we have an icon picker open?
                if($('.fa_icon_picker:visible').length){
                    $('.param_fa_icon').popover('hide');
                    if(prev_focus){
                        prev_focus.focus();
                    }
                }
            }
        });
    });

    // Parse initial JSON Schema
    try {
        schema = JSON.parse($('#json_schema').text());
    } catch(e){
        alert("Error - Schema JSON could not be parsed. See the browser console for details.");
        console.log(e);
    }

    // Build the schema builder
    try {
        $('#schema-builder').html( generate_obj(schema['properties'], 1) );
        init_group_sortable();
    } catch(e){
        alert("Error: Could not load schema JSON! Substituted for a blank schema.");
        console.error(e);
        schema = {
            "$schema": "http://json-schema.org/draft-07/schema",
            "$id": "https://raw.githubusercontent.com/YOUR_PIPELINE/master/nextflow_schema.json",
            "title": "Nextflow pipeline parameters",
            "description": "This pipeline uses Nextflow and processes some kind of data. The JSON Schema was built using the nf-core pipeline schema builder.",
            "type": "object",
            "properties": { }
        };
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    }

    // Add parameter button
    $('.add-param-btn').click(function(e){
        var new_id = 'new_param_'+new_param_idx;
        while (Object.keys(schema['properties']).indexOf(new_id) != -1) {
            new_param_idx += 1;
            new_id = 'new_param_'+new_param_idx;
        }
        var new_param = {
            "type": "string",
            "description": "",
            "default": ""
        };
        schema['properties'][new_id] = new_param;
        param_row = $( generate_param_row(new_id, new_param) );
        param_row.prependTo('#schema-builder').find('.param_id').select();
        scroll_to( param_row );
        schema_order_change();
        new_param_idx += 1;

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    // Add group button
    $('.add-group-btn').click(function(e){
        var new_id = 'New Group '+new_group_idx;
        while (Object.keys(schema['properties']).indexOf(new_id) != -1) {
            new_group_idx += 1;
            new_id = 'New Group '+new_group_idx;
        }
        var new_param = {
            "type": "object",
            "description": "",
            "default": "",
            "properties": {}
        };
        schema['properties'][new_id] = new_param;
        param_row = $( generate_group_row(new_id, new_param) );
        param_row.prependTo('#schema-builder').find('.param_id').select();
        scroll_to( param_row );
        init_group_sortable();
        schema_order_change();
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
                    new_schema['properties'] = {};
                    for(k in schema['properties']){
                        var new_k = k;
                        if(k == id){ new_k =  new_id};
                        new_schema['properties'][new_k] = schema['properties'][k];
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
            if($(this).val().trim().length == 0){
                delete new_param[param_key];
            } else {
                new_param[param_key] = $(this).val().trim();
            }

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
                        if(param_default.toLowerCase() == 'true'){
                            param['default'] = true;
                        } else {
                            param['default'] = false;
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

                // Convert bool default strings to bool
                if(param['type'] == 'boolean'){
                    if(!param.hasOwnProperty('default')){
                        param['default'] = false;
                    } else if(typeof param['default'] != "boolean"){
                        if(param['default'].toLowerCase() == 'true'){
                            param['default'] = true;
                        } else {
                            param['default'] = false;
                        }
                    }
                }

                // Convert number types
                if(['number', 'int', 'range'].indexOf(param['type']) != -1){
                    if(param.hasOwnProperty('default') && typeof param['default'] == "string" && param['default'].trim() != ''){
                        if(param['type'] == 'int'){
                            param['default'] = parseInt(param['default']);
                        } else {
                            param['default'] = parseFloat(param['default']);
                        }
                    }
                }

                update_param_in_schema(id, param);
            }

        }

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Keypress listener - move around with keyboard shortcuts
    //
    $('body').on('keydown', 'input, select', function(e){

        // Enter key
        if(e.which == 13){
            e.preventDefault();
            var current_class = '.'+$(this).attr('class').split(' ').join('.');
            var row = $(this).closest('.schema_row');
            // Shift + Enter - go up
            if(e.shiftKey){
                row.prev('.schema_row').find(current_class).focus();
            }
            // Just enter -go down
            else {
                row.next('.schema_row').find(current_class).focus();
            }
        }

        // Tab works as we want already by default

        // cmd + shift + , (open settings)
        if(e.which == 188 && e.shiftKey && (e.ctrlKey || e.metaKey)){
            // Get row of focussed element
            var row = $(':focus').closest('.schema_row');
            if(row.length == 1){
                var id = row.data('id');
                prev_focus = $(':focus');
                launch_settings_modal(id);
            }
        }

        // cmd + shif + up cursor (move row up)
        if(e.which == 38 && e.shiftKey && (e.ctrlKey || e.metaKey)){
            // Get row of focussed element
            var row = $(':focus').closest('.schema_row');
            if(row.hasClass('schema_group_row')){
                row = $(':focus').closest('.schema_group');
            }
            if(row.length == 1){
                prev_focus = $(':focus');
                var row_before = row.prev();
                row.insertBefore(row_before);
                schema_order_change();
                prev_focus.focus();
            }
        }

        // cmd + shif + down cursor (move row down)
        if(e.which == 40 && e.shiftKey && (e.ctrlKey || e.metaKey)){
            // Get row of focussed element
            var row = $(':focus').closest('.schema_row');
            if(row.hasClass('schema_group_row')){
                row = $(':focus').closest('.schema_group');
            }
            if(row.length == 1){
                prev_focus = $(':focus');
                var row_after = row.next();
                row.insertAfter(row_after);
                schema_order_change();
                prev_focus.focus();
            }
        }

        // Escape - hide icon picker
        if(e.which == 27){
            if($('.popover:visible').length){
                $('.param_fa_icon').popover('hide');
                prev_focus.focus();
            }
        }
    });

    //
    // Sorting - element has been moved
    //
    $('#schema-builder, .schema_group .card-body').on('sortstop', function(e, ui){
        schema_order_change();
    });

    function schema_order_change(){
        // Don't actually need to know where it landed - just rebuild schema from the DOM
        var new_schema = JSON.parse(JSON.stringify(schema));
        new_schema['properties'] = {};
        new_schema['required'] = [];
        $('.schema_row').each(function(idx, row){
            var id = $(row).data('id');
            var param = JSON.parse(JSON.stringify(find_param_in_schema(id)));

            // Check if we are inside a group
            if ($(this).parent('.card-body').length) {
                var group_id = $(this).parent().data('id');
                new_schema['properties'][group_id]['properties'][id] = param;
                if($(this).find('.param_required').is(':checked')){
                    new_schema['properties'][group_id]['required'].push(id);
                }
            } else {
                new_schema['properties'][id] = param;
                if($(this).find('.param_required').is(':checked')){
                    new_schema['required'].push(id);
                }
                // If a group, delete contents of that as well
                if(new_schema['properties'][id].hasOwnProperty('properties')){
                    new_schema['properties'][id]['properties'] = {};
                    new_schema['properties'][id]['required'] = [];
                }
            }
        });
        // Clear empty required arrays
        if(new_schema['required'].length == 0){
            delete new_schema['required'];
        }
        for(k in new_schema['properties']){
            if(new_schema['properties'][k].hasOwnProperty('required')){
                if(new_schema['properties'][k]['required'].length == 0){
                    delete new_schema['properties'][k]['required'];
                }
            }
        }

        schema = new_schema;
        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    };

    //
    // Required - required checkbox pressed
    //
    $('#schema-builder').on('change', 'input.param_required', function(){
        var id = $(this).closest('.schema_row').data('id');
        var is_required = $(this).is(':checked');
        set_required(id, is_required);
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
    // Help-text modal
    //
    $('#schema-builder').on('click', '.schema_row_help_text_icon', function(e){

        // Populate the help text modal
        var id = $(this).closest('.schema_row').data('id');
        var param = find_param_in_schema(id);
        $('#help_text_modal .modal-title').html('params.<span>'+id+'</span>');
        $('.helptext-preview-title').text('--'+id);
        if(param.description == undefined){
            param.description = '';
        }
        $('.helptext-preview-description').text(param.description);
        if(param.help_text == undefined){
            param.help_text = '';
        }
        $('#help_text_input').val(param.help_text);
        $('.helptext-preview-helptext').text(param.help_text);

        // Reset the tabs
        $('#help_text_modal .card-header-tabs .nav-link').removeClass('active');
        $('#help_text_modal .card-header-tabs .nav-link[href="#tab-helptext"]').addClass('active');
        $('#help_text_modal .tab-pane').removeClass('active show');
        $('#help_text_modal #tab-helptext').addClass('active show');

        // Show the modal
        prev_focus = $(this);
        $('#help_text_modal').modal('show');
    });

    // Focus the help text inputs
    $('#help_text_modal').on('shown.bs.modal', function(){
        $('#help_text_input').focus();
    });

    // Refocus page when modal closes
    $('#help_text_modal').on('hidden.bs.modal', function(){
        if(prev_focus){
            prev_focus.focus();
        }
    });

    // Re-render preview on tab change
    $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
        if($(e.target).attr('href') == '#tab-helptext-preview'){
            // Update basic text
            $('.helptext-preview-helptext').text($('#help_text_input').val());

            // Convert markdown to HTML
            var md_converter = new showdown.Converter();
            var help_text_html = md_converter.makeHtml( $('#help_text_input').val() );
            $('.helptext-html-preview .helptext-preview-helptext').html(help_text_html);
        }
    })

    // Save the help text
    $('#help_text_save').click(function(){
        var id = $('#help_text_modal .modal-title span').text();
        var param = find_param_in_schema(id);
        var help_text = $('#help_text_input').val();

        // Update the help-text icon
        var help_text_icon = help_text_icon_template;
        if(help_text.length > 0){
            help_text_icon = no_help_text_icon;
        }
        $(".schema_row[data-id='"+id+"'] .schema_row_help_text_icon i").replaceWith($(help_text_icon));

        // Update the schema
        if(help_text.length > 0){
            param.help_text = help_text;
        } else {
            delete param.help_text;
        }

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
        launch_settings_modal(id);
    });

    function launch_settings_modal(id){
        var param = find_param_in_schema(id);

        // Build modal
        var modal_header = 'params.<span>'+id+'</span>';
        if(param['type'] == 'object'){
            modal_header = '<span>'+id+'</span>';
        }
        $('#settings_modal .modal-title').html(modal_header);
        $('#settings_enum, #settings_pattern, #settings_minimum, #settings_maximum').val('');
        $('.settings_nothing_special, .settings_enum_group, .settings_pattern_group, .settings_minmax_group').hide();

        if(['boolean', 'object'].indexOf(param['type']) != -1){
            $('.settings_nothing_special').show();
        } else {
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

        $('#settings_modal').modal('show');
    }

    //
    // Settings Modal - save button
    //
    $('#settings_save').click(function(e){

        var id = $('#settings_modal .modal-title span').text();
        var param = find_param_in_schema(id);

        var settings = {};
        settings.pattern = $('#settings_pattern').val().trim();
        settings.minimum = $('#settings_minimum').val().trim();
        settings.maximum = $('#settings_maximum').val().trim();
        settings.enum = $('#settings_enum').val().trim().split('|');
        // Trim whitespace from each element and remove empties
        $.map(settings.enum, $.trim);
        settings.enum = settings.enum.filter(function (el) { return el.length > 0; });

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
            } else if(prev_focus){
                prev_focus.focus();
                prev_focus = false;
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

        for(k in schema['properties']){
            // Check if group
            if(schema['properties'][k].hasOwnProperty('properties')){
                for(j in schema['properties'][k]['properties']){
                    if(j == id){
                        delete schema['properties'][k]['properties'][j];
                    }
                    // Move contents of the group out if we're going to delete it
                    if(k == id){
                        // Move the HTML row out of the group
                        $('.schema_row[data-id="'+j+'"]').insertBefore(group_el);

                        // Move the schema param in to the top-level schema object
                        schema['properties'][j] = schema['properties'][k]['properties'][j];

                        // If it is required, set this on the top-level schema object
                        if(schema['properties'][k].hasOwnProperty('required') && schema['properties'][k]['required'].indexOf(j) != -1){
                            set_required(j, true);
                        }
                    }
                }
            }
            if(k == id){
                delete schema['properties'][k];
            }
        }

        row_el.remove();
        group_el.remove();

        // Update printed schema in page
        $('#json_schema').text(JSON.stringify(schema, null, 4));
    });

    //
    // Collapse group button
    //
    $('#schema-builder').on('click', '.schema_group_collapse', function(){
        $(this).closest('.schema_group').find('.card-body').slideToggle('fast');
        $(this).find('i').toggleClass('fa-angle-double-down fa-angle-double-up')
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

    var description = '';
    if(param['description'] != undefined){
        description = param['description'];
    }

    var default_input = '';
    if(param['type'] == 'boolean'){
        default_input = `
            <select class="param_key param_default" data-param_key="default">
                <option `+(param['default'] ? 'selected="selected"' : '')+`>True</option>
                <option `+(!param['default'] ? 'selected="selected"' : '')+`>False</option>
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
        default_input = '<div class="w-100"><input '+attrs+' class="param_key param_default" data-param_key="default"></div>';
    }

    var is_required = false;
    // Check that the required array exists
    if(schema.hasOwnProperty('required')){
        if (schema['required'].indexOf(id) !== -1) {
            is_required = true;
        }
    }
    // Lazy, just checking if it's in any group rather than specifically its own
    for(k in schema['properties']){
        if(schema['properties'][k].hasOwnProperty('required')){
            if (schema['properties'][k]['required'].indexOf(id) !== -1) {
                is_required = true;
            }
        }
    }

    var is_hidden = false;
    if(param['hidden']){
        is_hidden = true;
    }

    var fa_icon = '<i class="far fa-question-circle fa-fw param_fa_icon_missing"></i>';
    if(param['fa_icon'] != undefined && param['fa_icon'].trim().length > 0){
        var re = new RegExp('^fa[a-z -]+$');
        if(!re.test(param['fa_icon'])){
            console.error("FontAwesome icon did not match the regex: /^fa[a-z -]+$/ ('"+param['fa_icon']+"') - removing from schema.");
            delete param['fa_icon'];
            update_param_in_schema(id, param);
            $('#json_schema').text(JSON.stringify(schema, null, 4));
        } else {
            fa_icon = '<i class="'+param['fa_icon']+' fa-fw"></i>';
        }
    }

    var help_text_icon = help_text_icon_template;
    if(param['help_text'] != undefined && param['help_text'].trim().length > 0){
        help_text_icon = no_help_text_icon;
    }


    var results = `
    <div class="row schema_row border" data-id="`+id+`">
        <div class="col-auto align-self-center schema_row_grabber border-right">
            <i class="fas fa-grip-vertical"></i>
        </div>
        <button class="col-auto align-self-center param_fa_icon ">`+fa_icon+`</button>
        <div class="col schema-id">
            <label>ID
                <input type="text" class="text-monospace param_id" value="`+id+`">
            </label>
        </div>
        <div class="d-sm-none w-100"></div>
        <div class="col">
            <label>Description
                <input type="text" class="param_key param_description" data-param_key="description" value="`+description+`">
            </label>
        </div>
        <button class="col-auto align-self-center schema_row_help_text_icon">`+help_text_icon+`</button>
        <div class="col-auto">
            <label>Type
                <select class="param_key param_type" data-param_key="type">
                    <option `+(param['type'] == 'string' ? 'selected="selected"' : '')+` value="string">string</option>
                    <option `+(param['type'] == 'number' ? 'selected="selected"' : '')+` value="number">number</option>
                    <option `+(param['type'] == 'integer' ? 'selected="selected"' : '')+` value="integer">integer</option>
                    <option `+(param['type'] == 'range' ? 'selected="selected"' : '')+` value="range">range</option>
                    <option `+(param['type'] == 'boolean' ? 'selected="selected"' : '')+` value="boolean">boolean</option>
                </select>
            </label>
        </div>
        <div class="d-sm-none w-100"></div>
        <div class="col">
            <label>Default `+default_input+`</label>
        </div>
        <div class="col-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>R<span class="d-none d-lg-inline">equired</span>
                <input type="checkbox" `+(is_required ? 'checked="checked"' : '')+`" class="param_required">
            </label>
            `)+`
        </div>
        <div class="col-auto">
            `+(param['type'] == 'object' ? '' : `
            <label>H<span class="d-none d-lg-inline">idden</span>
                <input type="checkbox" `+(is_hidden ? 'checked="checked"' : '')+`" class="param_hidden">
            </label>
            `)+`
        </div>
        <div class="col-auto align-self-center schema_row_config border-left">
            <i class="fas fa-cog"></i>
        </div>
    </div>`;

    return results;
}

function generate_group_row(id, param, child_params){

    var description = '';
    if(param['description'] != undefined){
        description = param['description'];
    }

    if(child_params == undefined){
        child_params = '';
    }

    var fa_icon = '<i class="far fa-question-circle fa-fw param_fa_icon_missing"></i>';
    if(param['fa_icon'] != undefined && param['fa_icon'].trim().length > 0){
        var re = new RegExp('^fa[a-z -]+$');
        if(!re.test(param['fa_icon'])){
            console.error("FontAwesome icon did not match the regex: /^fa[a-z -]+$/ ('"+param['fa_icon']+"') - removing from schema.");
            delete param['fa_icon'];
            update_param_in_schema(id, param);
            $('#json_schema').text(JSON.stringify(schema, null, 4));
        } else {
            fa_icon = '<i class="'+param['fa_icon']+' fa-fw"></i>';
        }
    }

    var help_text_icon = help_text_icon_template;
    if(param['help_text'] != undefined && param['help_text'].trim().length > 0){
        help_text_icon = no_help_text_icon;
    }

    var results = `
    <div class="card schema_group" data-id="`+id+`">
        <div class="card-header p-0">
            <div class="row schema_row schema_group_row mb-0" data-id="`+id+`">
                <div class="col-auto align-self-center schema_row_grabber border-right">
                    <i class="fas fa-grip-vertical"></i>
                </div>
                <button class="col-auto align-self-center param_fa_icon ">`+fa_icon+`</button>
                <div class="col schema-id">
                    <label>Title
                        <input type="text" class="text-monospace param_id" value="`+id+`">
                    </label>
                </div>
                <button class="col-auto align-self-center schema_row_help_text_icon">`+help_text_icon+`</button>
                <div class="col">
                    <label>Description
                        <input type="text" class="param_key" data-param_key="description" value="`+description+`">
                    </label>
                </div>
                <div class="col-auto align-self-center schema_row_config border-left">
                    <i class="fas fa-cog"></i>
                </div>
                <div class="col-auto align-self-center schema_group_collapse">
                    <i class="fas fa-angle-double-down"></i>
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
        placeholder: 'schema_row_move_placeholder alert alert-secondary',
        connectWith: '.schema_group .card-body'
    });

    // Object Groups
    $(".schema_group .card-body").sortable({
        handle: '.schema_row_grabber',
        tolerance: 'pointer',
        placeholder: 'schema_row_move_placeholder alert alert-secondary',
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
    if(schema['properties'].hasOwnProperty(id)){
        num_hits += 1;
    }
    // Iterate through groups, looking for ID
    for(k in schema['properties']){
        if(schema['properties'][k].hasOwnProperty('properties')){
            if(schema['properties'][k]['properties'].hasOwnProperty(id)){
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
    if(!param.hasOwnProperty('default') || param['default'].length == 0){
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


function set_required(id, is_required){
    // Function to set the required flag in the JSON Schema for a given ID
    var schema_parent = null;
    // Get the parent object in the schema
    if(schema['properties'].hasOwnProperty(id)){
        schema_parent = schema;
    } else {
        // Iterate through groups, looking for ID
        for(k in schema['properties']){
            // Check if group
            if(schema['properties'][k].hasOwnProperty('properties')){
                if(schema['properties'][k]['properties'].hasOwnProperty(id)){
                    schema_parent =  schema['properties'][k];
                }
            }
        }
    }
    // Check that the required array exists
    if(!schema_parent.hasOwnProperty('required')){
        schema_parent['required'] = [];
    }
    if(is_required){
        schema_parent['required'].push(id);
    } else {
        var idx = schema_parent['required'].indexOf(id);
        if (idx !== -1) {
            schema_parent['required'].splice(idx, 1);
        }
    }
    // Remove required array if empty
    if(schema_parent['required'].length == 0){
        delete schema_parent['required'];
    }
    // Update printed schema in page
    $('#json_schema').text(JSON.stringify(schema, null, 4));
}


function find_param_in_schema(id){
    // Given an ID, find the param schema even if it's in a group
    // Assumes max one level of nesting and unique IDs everywhere

    // Simple case - not in a group
    if(schema['properties'].hasOwnProperty(id)){
        return schema['properties'][id];
    }

    // Iterate through groups, looking for ID
    for(k in schema['properties']){
        // Check if group
        if(schema['properties'][k].hasOwnProperty('properties')){
            if(schema['properties'][k]['properties'].hasOwnProperty(id)){
                return schema['properties'][k]['properties'][id];
            }
        }
    }

    console.warn("Could not find param '"+id+"'");
}

function update_param_in_schema(id, new_param){
    // Given an ID, find the param schema even if it's in a group
    // Assumes max one level of nesting and unique IDs everywhere

    // Simple case - not in a group
    if(schema['properties'].hasOwnProperty(id)){
        schema['properties'][id] = new_param;
        return true;
    }

    // Iterate through groups, looking for ID
    for(k in schema['properties']){
        // Check if group
        if(schema['properties'][k].hasOwnProperty('properties')){
            if(schema['properties'][k]['properties'].hasOwnProperty(id)){
                schema['properties'][k]['properties'][id] = new_param;
                return true;
            }
        }
    }

    console.warn("Could not find param '"+id+"'");
}
