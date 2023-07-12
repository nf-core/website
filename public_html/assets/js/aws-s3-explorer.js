// adapted from https://github.com/awslabs/aws-js-s3-explorer
// Copyright 2014-2018 Amazon.com, Inc. or its affiliates. All Rights Reserved.
// Licensed under the Apache License, Version 2.0 (the "License").
// You may not use this file except in compliance with the License. A copy
// of the License is located at
// https://aws.amazon.com/apache2.0/
// or in the "license" file accompanying this file. This file is distributed
// on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
// either express or implied. See the License for the specific language governing
// permissions and limitations under the License.

$(function () {
  AWS.config.region = s3exp_config.Region;
  $('#file-preview').hide();
  // Initialize S3 SDK and the moment library (for time formatting utilities)
  var ep = new AWS.Endpoint('cog.sanger.ac.uk');
  var s3 = new AWS.S3({endpoint: ep});
  moment().format();

  function bytesToSize(bytes) {
    var sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    if (bytes === 0) return '0 Bytes';
    var ii = parseInt(Math.floor(Math.log(bytes) / Math.log(1024)));
    return Math.round(bytes / Math.pow(1024, ii), 2) + ' ' + sizes[ii];
  }

  // Custom startsWith function for String prototype
  if (typeof String.prototype.startsWith != 'function') {
    String.prototype.startsWith = function (str) {
      return this.indexOf(str) == 0;
    };
  }

  // Custom endsWith function for String prototype
  if (typeof String.prototype.endsWith != 'function') {
    String.prototype.endsWith = function (str) {
      return this.slice(-str.length) == str;
    };
  }

  function object2hrefvirt(bucket, key) {
    var enckey = key
      .split('/')
      .map(function (x) {
        return encodeURIComponent(x);
      })
      .join('/');

    if (AWS.config.region === 'us-east-1') {
      return document.location.protocol + '//' + bucket + '.cog.sanger.ac.uk/' + enckey;
    } else {
      return document.location.protocol + '//' + bucket + '.cog.sanger.ac.uk/' + enckey;
    }
  }

  function object2hrefpath(bucket, key) {
    var enckey = key
      .split('/')
      .map(function (x) {
        return encodeURIComponent(x);
      })
      .join('/');

    if (AWS.config.region === 'us-east-1') {
      return document.location.protocol + '//.cog.sanger.ac.uk/' + bucket + '/' + enckey;
    } else {
      return document.location.protocol + '.cog.sanger.ac.uk/' + bucket + '/' + enckey;
    }
  }

  function sanitize_html(string) {
    const map = {
      '&': '&amp;',
      '<': '&lt;',
      '>': '&gt;',
      '"': '&quot;',
      "'": '&#x27;',
      '`': '&grave;',
      '/': '&#x2F;',
    };
    const reg = /[&<>"'`/]/gi;
    return string.replace(reg, (match) => map[match]);
  }

  function isthisdocument(bucket, key) {
    return key === 'index.html';
  }

  function isfolder(path) {
    return path.endsWith('/');
  }

  // Convert cars/vw/golf.png to golf.png
  function fullpath2filename(path) {
    return sanitize_html(path.replace(/^.*[\\\/]/, ''));
  }

  // Convert cars/vw/golf.png to cars/vw
  function fullpath2pathname(path) {
    return sanitize_html(path.substring(0, path.lastIndexOf('/')));
  }

  // Convert cars/vw/ to vw/
  function prefix2folder(prefix) {
    var parts = prefix.split('/');
    return sanitize_html(parts[parts.length - 2] + '/');
  }

  // Remove hash from document URL
  function removeHash() {
    history.pushState('', document.title, window.location.pathname + window.location.search);
  }
  function fetch_file_size(url) {
    var file_size = fetch(url, { method: 'HEAD' }).then(function (resp) {
      return resp.headers.get('content-length');
    });
    return file_size;
  }
  function fetch_preview(url, file_size) {
    var extension = url.substring(url.lastIndexOf('.'), url.length);
    var split_url = url.split('/');
    var filename = split_url.slice(-1)[0];
    var s3_url = 's3://' + split_url[2].split('.')[0] + '/' + split_url.slice(3).join('/');
    var btn_group = `<div class="btn-group" role="group">
                                <a class="btn btn-outline-secondary download-file-btn" href="${url}" target="_blank">Download file</a>
                                <button type="button" class="btn btn-outline-secondary copy-url" data-target=${url}>Copy URL</button>
                                <button type="button" class="btn btn-outline-secondary copy-url" data-target=${s3_url}>Copy S3 URL</button>
                            </div>`;
    var header =
      '<div class="card-header d-flex justify-content-between align-items-center">' + filename + btn_group + '</div>';
    if (file_size > 10000000) {
      data =
        '<div class="alert alert-warning text-center mb-0" role="alert"><i class="fad fa-exclamation-triangle"></i> The file is too big to be previewed.</div>';
      $('#file-preview').show();
      $('#file-preview').html(header + '<div class="card-body">' + data + '</div>');
      var el_offset = $('#file-preview').offset().top - 140;
      $([document.documentElement, document.body]).animate({ scrollTop: el_offset }, 500);
      return;
    }

    fetch(url)
      .then(function (response) {
        if (response.status !== 200) {
          console.log('Looks like there was a problem. Status Code: ' + response.status);
          return;
        }
        // Examine the text in the response
        response.text().then(function (data) {
          if (['.bam', '.bai', '.gz', '.zip'].includes(extension)) {
            data =
              '<div class="alert alert-warning text-center mb-0" role="alert"><i class="fad fa-exclamation-triangle"></i> No preview available for binary or compressed files.</div>';
          } else if (!['.html', '.pdf', '.png', '.jpg', '.jpeg', '.svg'].includes(extension)) {
            var lang = ['.json', '.yaml', '.yml'].includes(extension) ? extension.slice(1) : 'plaintext';
            data = '<pre><code class="language-' + lang + '">' + sanitize_html(data) + '</code></pre>';
          } else if (extension === '.html') {
            data =
              '<iframe srcdoc="' + sanitize_html(data) + '" style="border:none; width:100%; height:1000px;"></iframe>';
          } else if (extension === '.pdf') {
            data = `<object data="${url}" type="application/pdf" style="border:none; width:100%; height:1000px;">
                                            <embed src="${url}" type="application/pdf" />
                                        </object>`;
          } else if (['.png', '.jpg', '.jpeg', '.svg'].includes(extension)) {
            data = '<img src="' + response.url + '"/>';
          }
          if (data !== '<pre><code class="hljs language-' + lang + '"></code></pre>') {
            $('#file-preview').show();
            $('#file-preview').html(header + '<div class="card-body">' + data + '</div>');
            hljs.highlightAll();
            var el_offset = $('#file-preview').offset().top - 140;
            $([document.documentElement, document.body]).animate({ scrollTop: el_offset }, 500);
          }
        });
      })
      .catch(function (err) {
        console.log('Fetch Error :-S', err);
      });
  }

  $('body').on('click', '.download-file-btn ', function () {
    if ($(this).attr('href') && $(this).attr('href').length > 0) {
      var url = $(this).attr('href');
    } else {
      var url = $(this).siblings('div').children('a')[0].href;
    }
    window.open(url, '_blank');
  });

  //copy text to clipboard
  $('.toast').toast();

  $('body').on('click', '.copy-url', function (e) {
    var text = e.currentTarget.dataset.target;
    var $tmp = $('<input>');

    $('body').append($tmp);
    $tmp.val(text).select();
    document.execCommand('copy');
    $tmp.remove();
    $('#url-copied').toast('show');
  });

  //update view when url-hash changes
  $(window).on('hashchange', function (e) {
    if (!prefix.endsWith('/')) {
      prefix = window.location.hash.substring(1);
      prefix = prefix.substring(0, prefix.lastIndexOf('/') + 1);
      if (window.location.hash.split('/').length > 2) {
        s3exp_config['Prefix'] = prefix;
      }
      (s3exp_lister = s3list(s3exp_config, s3draw)).go();
    }
  });

  function folder2breadcrumbs(data) {
    // console.log('Bucket: ' + data.params.Bucket);
    // console.log('Prefix: ' + data.params.Prefix);

    if (data.params.Prefix && data.params.Prefix.length > 0) {
      // console.log('Set hash: ' + data.params.Prefix)
      window.location.hash = data.params.Prefix;
      if (s3exp_config.Suffix !== undefined && s3exp_config.Suffix !== '' && window.location.hash.endsWith('/')) {
        window.location.hash += s3exp_config.Suffix;
      }
    } else {
      // console.log('Remove hash');
      removeHash();
    }

    // The parts array will contain the bucket name followed by all the
    // segments of the prefix, exploded out as separate strings.
    var parts = [data.params.Bucket];

    if (data.params.Prefix) {
      parts.push.apply(
        parts,
        data.params.Prefix.endsWith('/') ? data.params.Prefix.slice(0, -1).split('/') : data.params.Prefix.split('/'),
      );
    }

    // console.log('Parts: ' + parts + ' (length=' + parts.length + ')');

    // Empty the current breadcrumb list
    $('#breadcrumb li').remove();

    // Now build the new breadcrumb list
    var buildprefix = '';
    $.each(parts, function (ii, part) {
      var ipart;

      // Add the bucket (the bucket is always first)
      if (ii === 0 || ii === 1) {
        var a1 = $('<span>').text(part);
        ipart = $('<li>').addClass('breadcrumb-item').append(a1);
        if (ii === 1) {
          buildprefix += part + '/';
        }
        // Else add the folders within the bucket
      } else {
        buildprefix += part + '/';

        if (ii == parts.length - 1) {
          ipart = $('<li>').addClass('breadcrumb-item text-break active').text(part);
        } else {
          var a2 = $('<a>').attr('href', '#').append(part);
          ipart = $('<li>').addClass('breadcrumb-item text-break').append(a2);

          // Closure needed to enclose the saved S3 prefix
          (function () {
            var saveprefix = buildprefix;
            // console.log('Part: ' + part + ' has buildprefix: ' + saveprefix);
            a2.click(function (e) {
              e.preventDefault();
              // console.log('Breadcrumb click object prefix: ' + saveprefix);
              s3exp_config = {
                Bucket: data.params.Bucket,
                Prefix: saveprefix,
                Delimiter: data.params.Delimiter,
              };
              (s3exp_lister = s3list(s3exp_config, s3draw)).go();
            });
          })();
        }
      }
      $('#breadcrumb').append(ipart);
    });
    $('.title-bar .copy-url')[0].dataset.target = 's3://' + parts.join('/');
  }

  function s3draw(data, complete) {
    $('li.li-bucket').remove();
    $('#file-preview').hide();

    folder2breadcrumbs(data);
    if (data.Contents.length > 0 || data.CommonPrefixes.length > 0) {
      var path = data.params.Prefix.split('/');
      if (path.length > 3) {
        $('#tb-s3objects')
          .DataTable()
          .rows.add([
            {
              Key: path.slice(0, path.length - 2).join('/') + '/',
              render_name: false,
            },
          ]);
      }
      // Add each part of current path (S3 bucket plus folder hierarchy) into the table
      $.each(data.CommonPrefixes, function (i, prefix) {
        $('#tb-s3objects')
          .DataTable()
          .rows.add([
            {
              Key: prefix.Prefix,
              render_name: true,
            },
          ]);
      });

      // Add S3 objects to DataTable
      $('#tb-s3objects').DataTable().rows.add(data.Contents).draw();
    } else {
      $('#tb-s3objects')
        .DataTable()
        .rows.add([
          {
            Key: '/',
            render_name: false,
          },
        ]);
      $('#tb-s3objects').DataTable().rows.add(data.Contents).draw();
    }
  }

  function s3list(config, completecb) {
    // console.log('s3list config: ' + JSON.stringify(config));
    var params = {
      Bucket: config.Bucket,
      Prefix: config.Prefix,
      Delimiter: config.Delimiter,
    };
    var scope = {
      Contents: [],
      CommonPrefixes: [],
      params: params,
      stop: false,
      completecb: completecb,
    };

    return {
      // This is the callback that the S3 API makes when an S3 listObjectsV2
      // request completes (successfully or in error). Note that a single call
      // to listObjectsV2 may not be enough to get all objects so we need to
      // check if the returned data is truncated and, if so, make additional
      // requests with a 'next marker' until we have all the objects.
      cb: function (err, data) {
        if (err) {
          console.log('Error: ' + JSON.stringify(err));
          console.log('Error: ' + err.stack);
          scope.stop = true;
          $('#bucket-loader').removeClass('fa-spin');
          alert('Error accessing S3 bucket ' + scope.params.Bucket + '. Error: ' + err);
        } else {
          // console.log('Data: ' + JSON.stringify(data));
          // console.log("Options: " + $("input[name='optionsdepth']:checked").val());

          // Store marker before filtering data
          if (data.IsTruncated) {
            if (data.NextContinuationToken) {
              scope.params.ContinuationToken = data.NextContinuationToken;
            }
          }
          // Filter the folders out of the listed S3 objects
          // (could probably be done more efficiently)
          // console.log("Filter: remove folders");
          data.Contents = data.Contents.filter(function (el) {
            return el.Key !== scope.params.Prefix;
          });

          // Optionally, filter the root index.html out of the listed S3 objects
          if (HIDE_INDEX) {
            // console.log("Filter: remove index.html");
            data.Contents = data.Contents.filter(function (el) {
              return el.Key !== 'index.html';
            });
          }

          // Accumulate the S3 objects and common prefixes
          scope.Contents.push.apply(scope.Contents, data.Contents);
          scope.CommonPrefixes.push.apply(scope.CommonPrefixes, data.CommonPrefixes);

          // Update badge count to show number of objects read
          // $('#badgecount').text(scope.Contents.length + scope.CommonPrefixes.length);

          if (scope.stop) {
            // console.log('Bucket ' + scope.params.Bucket + ' stopped');
          } else if (data.IsTruncated) {
            // console.log('Bucket ' + scope.params.Bucket + ' truncated');
            s3.makeUnauthenticatedRequest('listObjectsV2', scope.params, scope.cb);
          } else {
            // console.log('Bucket ' + scope.params.Bucket + ' has ' + scope.Contents.length + ' objects, including ' + scope.CommonPrefixes.length + ' prefixes');
            delete scope.params.ContinuationToken;
            if (scope.completecb) {
              scope.completecb(scope, true);
            }
            $('#bucket-loader').removeClass('fa-spin');
          }
        }
      },

      // Start the spinner, clear the table, make an S3 listObjectsV2 request
      go: function () {
        scope.cb = this.cb;
        $('#bucket-loader').addClass('fa-spin');
        $('#tb-s3objects').DataTable().clear();
        s3.makeUnauthenticatedRequest('listObjectsV2', scope.params, this.cb);
      },

      stop: function () {
        scope.stop = true;
        delete scope.params.ContinuationToken;
        if (scope.completecb) {
          scope.completecb(scope, false);
        }
        $('#bucket-loader').removeClass('fa-spin');
      },
    };
  }

  function resetDepth() {
    $('#tb-s3objects').DataTable().column(1).visible(false);
    $('input[name="optionsdepth"]').val(['folder']);
    $('input[name="optionsdepth"][value="bucket"]').parent().removeClass('active');
    $('input[name="optionsdepth"][value="folder"]').parent().addClass('active');
  }

  $(document).ready(function () {
    // console.log('ready');

    function renderObject(data, type, full) {
      if (isthisdocument(s3exp_config.Bucket, data)) {
        // console.log("is this document: " + data);
        return fullpath2filename(data);
      } else if (isfolder(data)) {
        // console.log("is folder: " + data);
        if (full.render_name) {
          return (
            '<a data-s3="folder" data-prefix="' +
            sanitize_html(data) +
            '" href="' +
            object2hrefvirt(s3exp_config.Bucket, data) +
            '"><i class="fas fa-folder"></i> ' +
            prefix2folder(data) +
            '</a>'
          );
        } else {
          if (data === '/') {
            data = s3exp_config.Prefix.split('/').slice(0, -2).join('/') + '/';
          }
          return (
            '<div class="d-flex justify-content-between align-items-center"><div><a data-s3="folder" data-prefix="' +
            sanitize_html(data) +
            '" href="' +
            object2hrefvirt(s3exp_config.Bucket, data) +
            '"><i class="fad fa-folder-open"></i> ..</a></div></div>'
          );
        }
      } else {
        var icon = '<i class="fas fa-file"></i> ';

        // console.log("not folder/this document: " + data);
        return (
          '<div class="d-flex justify-content-between align-items-center"><div><a data-s3="object" href="' +
          object2hrefvirt(s3exp_config.Bucket, data) +
          '"download="' +
          fullpath2filename(data) +
          '" data-size=' +
          full.Size +
          '>' +
          icon +
          fullpath2filename(data) +
          '</a></div></div>'
        );
      }
    }

    function renderFolder(data, type, full) {
      return isfolder(data) ? '' : fullpath2pathname(data);
    }

    // Initial DataTable settings
    $('#tb-s3objects').DataTable({
      info: false,
      paging: false,
      searching: false,
      iDisplayLength: 50,
      language: { emptyTable: 'No data available' },
      order: [
        [1, 'asc'],
        [0, 'asc'],
      ],
      aoColumnDefs: [
        {
          aTargets: [0],
          mData: 'Key',
          mRender: function (data, type, full) {
            return type == 'display' ? renderObject(data, type, full) : data;
          },
          sType: 'key',
        },
        {
          aTargets: [1],
          mData: 'Key',
          width: '68%',
          mRender: function (data, type, full) {
            return renderFolder(data, type, full);
          },
        },
        {
          aTargets: [2],
          mData: 'LastModified',
          width: '20%',
          className: 'text-end',
          mRender: function (data, type, full) {
            return data ? moment(data).fromNow() : '';
          },
        },
        {
          aTargets: [3],
          width: '12%',
          className: 'text-end',
          mData: function (data, type, val) {
            return data.Size ? (type == 'display' ? bytesToSize(data.Size) : data.Size) : '';
          },
        },
      ],
    });
    $('#tb-s3objects').DataTable().column(s3exp_columns.key).visible(false);
    $('#tb-s3objects').wrap('<div class="table-responsive-md"></div>'); // Make table responsive

    // Custom sort for the Key column so that folders appear before objects
    $.fn.dataTableExt.oSort['key-asc'] = function (a, b) {
      var x = (isfolder(a) ? '0-' + a : '1-' + a).toLowerCase();
      var y = (isfolder(b) ? '0-' + b : '1-' + b).toLowerCase();
      return x < y ? -1 : x > y ? 1 : 0;
    };

    $.fn.dataTableExt.oSort['key-desc'] = function (a, b) {
      var x = (isfolder(a) ? '1-' + a : '0-' + a).toLowerCase();
      var y = (isfolder(b) ? '1-' + b : '0-' + b).toLowerCase();
      return x < y ? 1 : x > y ? -1 : 0;
    };

    // Delegated event handler for S3 object/folder clicks. This is delegated
    // because the object/folder rows are added dynamically and we do not want
    // to have to assign click handlers to each and every row.
    $('#tb-s3objects').on('click', 'a', function (event) {
      event.preventDefault();
      var target = event.target.href ? event.target : $(event.target).parent('a')[0];
      let hash = window.location.hash;
      hash = hash.substring(1, hash.lastIndexOf('/') + 1) + target.download;
      window.location.hash = hash;

      // console.log("target href=" + target.href);
      // console.log("target dataset=" + JSON.stringify(target.dataset));

      // If the user has clicked on a folder then navigate into that folder
      if (target.dataset.s3 === 'folder') {
        resetDepth();
        delete s3exp_config.ContinuationToken;
        s3exp_config.Prefix = target.dataset.prefix;
        s3exp_config.Delimiter = '/';
        (s3exp_lister = s3list(s3exp_config, s3draw)).go();
        // Else user has clicked on an object so preview it underneath
      } else {
        let url = target.href;
        let file_size = target.dataset.size;
        $('tr.table-active').removeClass('table-active');
        $(this).parents('tr').addClass('table-active');
        fetch_preview(url, file_size);
      }
      return false;
    });

    if (s3exp_config.Bucket) {
      if (s3exp_config.Prefix.endsWith('results-/')) {
        s3exp_config.Prefix = s3exp_config.Prefix.replace('results-/', '');
      }
      (s3exp_lister = s3list(s3exp_config, s3draw)).go(); //initial load
    }
  });
});
