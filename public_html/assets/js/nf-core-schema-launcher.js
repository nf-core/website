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

    // Parse initial JSON Schema
    try {
        schema = JSON.parse($('#json_schema').text());
    } catch(e){
        alert("Error - Schema JSON could not be parsed. See the browser console for details.");
        console.log(e);
    }

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
