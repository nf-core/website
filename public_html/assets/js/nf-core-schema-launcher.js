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
    var validation_error = false;
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
        if(!validation_error){
            $('#progress_section').html(newLabel);
        }
    };

    // Page-scroll links
    ///////////////////////////////////////////////////
    // TODO - NOT WORKING
    ///////////////////////////////////////////////////
    $('body').on('.scroll_to_link', 'click', function(e){
        e.preventDefault();
        scroll_to($(this));
    });

    // Validate form on submit - snippet from Bootstrap docs
    // Fetch all the forms we want to apply custom Bootstrap validation styles to
    var forms = document.getElementsByClassName('needs-validation');
    $('#form_validation_error_toast').toast({ 'autohide': false });
    // Loop over them and prevent submission
    var validation = Array.prototype.filter.call(forms, function(form) {
      form.addEventListener('submit', function(event) {
        if (form.checkValidity() === false) {
            event.preventDefault();
            event.stopPropagation();
            validation_error = true;
            $('.btn-launch').removeClass('btn-primary').addClass('btn-danger');
            $('.progress-bar').addClass('bg-danger');
            $('.validation-warning').show();
            $('#form_validation_error_toast').toast('show');
            $('#progress_section').html("Error validating inputs");
            $('#validation_fail_list').html('');
            $('input:invalid').each(function(){
                $('#validation_fail_list').append('<li><a href="#'+$(this).attr('id')+'" class="scroll_to_link"><code>'+$(this).attr('name').replace('params_', '--').replace('nxf_flag_', '')+'</code></a></li>')
            })
            scroll_to($('input:invalid').first());
        }
        form.classList.add('was-validated');

        // Radio button form group classes
        $('.form-control:has(.form-check-input[type="radio"]:valid)').addClass('radio-form-control-valid');
        $('.form-control:has(.form-check-input[type="radio"]:invalid)').addClass('radio-form-control-invalid');

        // Input group addons
        $('.input-group:has(input:valid) .input-group-prepend, .input-group:has(input:valid) .input-group-append').addClass('input-group-append-valid');
        $('.input-group:has(input:invalid) .input-group-prepend, .input-group:has(input:invalid) .input-group-append').addClass('input-group-append-invalid');

      }, false);
    });

});
