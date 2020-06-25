//
// nf-core-schema-launcher.js
// Custom javascript for the nf-core JSON Schema Launcher interface
//

$(function () {

    // Show the cache expiration time in local timezone
    $('.cache_expires_at span').text(
        moment.unix($('.cache_expires_at span').text()).calendar().replace(/^\w/, function (chr) { return chr.toLowerCase(); })
    );
    $('.cache_expires_at').show();

    // Show / hide hidden fields
    $('.btn-show-hidden-fields').click(function(){
        $('.is_hidden, .is_not_hidden').toggleClass('is_hidden is_not_hidden');
    });

    // Listen to the page scroll
    window.onscroll = function(){
        // Progress bar
        var winScroll = document.body.scrollTop || document.documentElement.scrollTop;
        var formTop = document.getElementById("schema_launcher_form").offsetTop;
        var height = document.documentElement.scrollHeight - document.documentElement.clientHeight - formTop;
        var scrolled = ((winScroll - formTop) / height) * 100;
        if(winScroll < formTop){ scrolled = 0; }
        $('.progress-bar').css('width', scrolled+"%").attr('area-valuenow', scrolled);

        // Jump to section dropdown
        var newLabel = 'Jump to section';
        $('legend:visible').each(function(){
            var this_offset = $(this).closest('fieldset').offset().top - 30;
            if(winScroll > this_offset && winScroll < this_offset + $(this).closest('fieldset').outerHeight(true)){
                newLabel = $(this).text();
            }
        });
        $('#dropdownMenuButton span').text(newLabel);
    };

    // Parse initial JSON Schema
    try {
        schema = JSON.parse($('#json_schema').text());
    } catch(e){
        alert("Error - Schema JSON could not be parsed. See the browser console for details.");
        console.log(e);
    }

    // Validate form on submit
    $('#schema_launcher_form').on('submit', function(event) {
        if (form.checkValidity() === false) {
            event.preventDefault();
            event.stopPropagation();
        }
        form.classList.add('was-validated');
    });

    //
    // FINISHED button
    //
    // Toggle between panels
    $('.launcher-panel-btn').click(function(){

        // Post the results to PHP when finished
        if($(this).data('target') == '#params-finished'){
            $('#schema-send-status').text("Saving schema..");

            post_data = {
                'post_content': 'json_schema_launcher',
                'version': 'web_builder',
                'status': 'launch_params_complete',
                'api': 'true',
                'cache_id': $('#params_cache_id').text(),
                'schema': JSON.stringify(schema),
                'nxf_flags': JSON.stringify({
                    "-resume": true,
                    "-revision": "dev"
                }),
                'input_params': JSON.stringify({
                    "input": "./design.csv",
                    "broad_cutoff": "0.2",
                    "fasta": "foobar"
                })
            };
            $.post( "json_schema_launch", post_data).done(function( returned_data ) {
                console.log("Sent schema to API. Response:", returned_data);
                if(returned_data.status == 'recieved'){
                    $('#schema-send-status').text("Ok, that's it - done!");
                } else {
                    console.log("Data sent to API:", post_data);
                    $('#schema-send-status').text("Oops, something went wrong!");
                }
            });
        }
    });



});
