//
// nf-core-schema-launcher.js
// Custom javascript for the nf-core JSON Schema Launcher interface
//

// Globals
validation_error = false;
section_label = 'Nextflow command-line flags';

$(function () {
  // Show the cache expiration time in local timezone
  $('.cache_expires_at span').text(
    moment
      .unix($('.cache_expires_at span').text())
      .calendar()
      .replace(/^\w/, function (chr) {
        return chr.toLowerCase();
      })
  );
  $('.cache_expires_at').show();

  // Show / hide hidden fields
  $('.btn-show-hidden-fields').click(function (e) {
    e.preventDefault();
    $('.is_hidden, .is_not_hidden').toggleClass('is_hidden is_not_hidden');
  });

  // Launch select-pipeline form
  $('#launch-pipeline-name').change(function () {
    var wf_name = $(this).val();
    var option_el = $(this).find(':selected');
    var releases = option_el.data('releases');
    if (wf_name == '') {
      $('#launch-pipeline-release').html('<option>Pipeline release</option>');
      $('#launch-pipeline-release').attr('disabled', true);
      $('#launch-pipeline-submit').attr('disabled', true);
    } else {
      $('#launch-pipeline-release').html('').attr('disabled', false);
      for (idx in releases) {
        $('#launch-pipeline-release').append('<option>' + releases[idx] + '</option>');
      }
      $('#launch-pipeline-submit').attr('disabled', false);
    }
  });

  // Listen to the page scroll
  if (document.getElementById('schema_launcher_form')) {
    window.onscroll = function () {
      // Progress bar width
      var winScroll = document.body.scrollTop || document.documentElement.scrollTop;
      var formTop = document.getElementById('schema_launcher_form').offsetTop;
      var height = document.documentElement.scrollHeight - document.documentElement.clientHeight - formTop;
      var scrolled = ((winScroll - formTop) / height) * 100;
      if (winScroll < formTop) {
        scrolled = 0;
      }

      // Jump to section dropdown
      $('legend:visible').each(function () {
        var this_offset = $(this).closest('fieldset').offset().top - 30;
        if (winScroll > this_offset && winScroll < this_offset + $(this).closest('fieldset').outerHeight(true)) {
          section_label = $(this).text();
        }
      });

      // Update progress bar
      $('.progress-bar')
        .css('width', scrolled + '%')
        .attr('area-valuenow', scrolled);
      if (!validation_error) {
        $('#progress_section').html(section_label);
      }
    };
  }

  // Page-scroll links
  $('body').on('click', '.scroll_to_link', function (e) {
    e.preventDefault();
    scroll_to($($(this).attr('href')), 140);
  });

  // Validate form on submit - snippet from Bootstrap docs
  $('#form_validation_error_toast').toast({ autohide: false });
  $('#schema_launcher_form').on('submit change keyup', function (e) {
    // Only run validation on change and keyup if we have already submitted / validated
    if (e.type == 'change' || e.type == 'keyup') {
      if (!$(this).hasClass('was-validated')) {
        return;
      }
    }

    // Validate form
    if ($(this)[0].checkValidity() === false) {
      e.preventDefault();
      e.stopPropagation();
      validation_error = true;
      set_validation_styles(false);
      // If form was submitted, scroll to first error
      if (e.type == 'submit') {
        scroll_to($('input:invalid, select:invalid').first(), 140);
        $('input:invalid, select:invalid').first().focus();
      }
    } else {
      validation_error = false;
      set_validation_styles(true);
    }
    $(this).addClass('was-validated');

    // Radio button form group classes
    $('.form-control:has(.form-check-input[type="radio"]:valid)').addClass('radio-form-control-valid');
    $('.form-control:has(.form-check-input[type="radio"]:invalid)').addClass('radio-form-control-invalid');

    // Input group addons
    $('.input-group:has(input:valid)').addClass('input-group-valid');
    $('.input-group:has(input:invalid, select:invalid)').addClass('input-group-invalid');
  });
});

function set_validation_styles(form_is_valid) {
  // Reset
  $('#validation_fail_list').html('');
  $('.form-control:has(.form-check-input[type="radio"])').removeClass(
    'radio-form-control-valid radio-form-control-invalid'
  );
  $('.input-group .input-group-prepend').removeClass('input-group-valid input-group-invalid');

  // Invalid
  if (!form_is_valid) {
    var invalid_els = $('input:invalid, select:invalid');
    $('.btn-launch').removeClass('btn-primary').addClass('btn-danger');
    $('.progress-bar').addClass('bg-danger');
    $('.validation-warning').show();
    $('#form_validation_error_toast').toast('show');
    $('#progress_section')
      .html(invalid_els.length + ' parameter' + (invalid_els.length > 1 ? 's' : '') + ' with errors')
      .removeClass('text-muted')
      .addClass('text-danger');
    invalid_els.each(function () {
      $('#validation_fail_list').append(
        '<li><a href="#' +
          $(this).attr('id') +
          '" class="scroll_to_link"><code>' +
          $(this).attr('name').replace('params_', '--').replace('nxf_flag_', '') +
          '</code></a></li>'
      );
    });
  }
  // Valid
  else {
    $('.btn-launch').addClass('btn-primary').removeClass('btn-danger');
    $('.progress-bar').removeClass('bg-danger');
    $('.validation-warning').hide();
    $('#form_validation_error_toast').toast('hide');
    $('#progress_section').html(section_label).addClass('text-muted').removeClass('text-danger');
  }
}
