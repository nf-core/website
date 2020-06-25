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
        // Progress bar width
        var winScroll = document.body.scrollTop || document.documentElement.scrollTop;
        var formTop = document.getElementById("schema_launcher_form").offsetTop;
        var height = document.documentElement.scrollHeight - document.documentElement.clientHeight - formTop;
        var scrolled = ((winScroll - formTop) / height) * 100;
        if(winScroll < formTop){ scrolled = 0; }

        // Jump to section dropdown
        var newLabel = 'Nextflow command-line flags';
        $('legend:visible').each(function(){
            var this_offset = $(this).closest('fieldset').offset().top - 30;
            if(winScroll > this_offset && winScroll < this_offset + $(this).closest('fieldset').outerHeight(true)){
                newLabel = $(this).text();
            }
        });

        // Update progress bar
        $('.progress-bar').css('width', scrolled+"%").attr('area-valuenow', scrolled);
        $('#progress_section').html(newLabel);
    };

    // Parse initial JSON Schema
    try {
        schema = JSON.parse($('#json_schema').text());
    } catch(e){
        alert("Error - Schema JSON could not be parsed. See the browser console for details.");
        console.log(e);
    }

    // Validate form on submit - snippet from Bootstrap docs
    // Fetch all the forms we want to apply custom Bootstrap validation styles to
    var forms = document.getElementsByClassName('needs-validation');
    // Loop over them and prevent submission
    var validation = Array.prototype.filter.call(forms, function(form) {
      form.addEventListener('submit', function(event) {
        if (form.checkValidity() === false) {
          event.preventDefault();
          event.stopPropagation();
        }
        form.classList.add('was-validated');
      }, false);
    });

});
