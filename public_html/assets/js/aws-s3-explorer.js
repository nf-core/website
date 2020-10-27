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
    $("#file-preview").hide();
    // Initialize S3 SDK and the moment library (for time formatting utilities)
    var s3 = new AWS.S3();
    moment().format();

    function bytesToSize(bytes) {
    var sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    if (bytes === 0) return '0 Bytes';
    var ii = parseInt(Math.floor(Math.log(bytes) / Math.log(1024)));
    return Math.round(bytes / Math.pow(1024, ii), 2) + ' ' + sizes[ii];
    }

    // Custom startsWith function for String prototype
    if (typeof String.prototype.startsWith != 'function') {
    String.prototype.startsWith = function(str) {
        return this.indexOf(str) == 0;
    };
    }

    // Custom endsWith function for String prototype
    if (typeof String.prototype.endsWith != 'function') {
    String.prototype.endsWith = function(str) {
        return this.slice(-str.length) == str;
    };
    }

    function object2hrefvirt(bucket, key) {
    var enckey = key.split('/').map(function(x) { return encodeURIComponent(x); }).join('/');

    if (AWS.config.region === "us-east-1") {
        return document.location.protocol + '//' + bucket + '.s3.amazonaws.com/' + enckey;
    } else {
        return document.location.protocol + '//' + bucket + '.s3-' + AWS.config.region + '.amazonaws.com/' + enckey;
    }
    }

    function object2hrefpath(bucket, key) {
    var enckey = key.split('/').map(function(x) { return encodeURIComponent(x); }).join('/');

    if (AWS.config.region === "us-east-1") {
        return document.location.protocol + "//s3.amazonaws.com/" + bucket + "/" + enckey;
    } else {
        return document.location.protocol + "//s3-" + AWS.config.region + ".amazonaws.com/" + bucket + "/" + enckey;
    }
    }

    function sanitize_html(string) {
            const map = {
                '&': '&amp;',
                '<': '&lt;',
                '>': '&gt;',
                '"': '&quot;',
                "'": '&#x27;',
                "`": '&grave;',
                "/": '&#x2F;',
            };
            const reg = /[&<>"'`/]/ig;
            return string.replace(reg, (match) => (map[match]));
        }

    function isthisdocument(bucket, key) {
    return key === "index.html";
    }

    function isfolder(path) {
        return path.endsWith('/');
    }
    function haspreview(data) {
        var extension = data.Key.substring(data.Key.lastIndexOf('.'),data.Key.length);
        var preview = [".csv",".tsv",".out",".table",".txt",".html"].includes(extension) && data.Size<80000000;
        return preview
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
        history.pushState("", document.title, window.location.pathname + window.location.search);
    }

    $("body").on("click",".preview-file-btn ",function(e){
        var url = $(this).siblings('div').children('a')[0].href;
        var extension = url.substring(url.lastIndexOf('.'),url.length);
        var split_url = url.split('/');
        var filename = split_url.slice(-1)[0];
        var s3_url = "s3://"+ split_url[2].split(".")[0] + "/" + split_url.slice(3).join("/");

        fetch(url)
            .then(
                function(response) {
                    if (response.status !== 200) {
                        console.log('Looks like there was a problem. Status Code: ' + response.status);
                        return;
                    }
                    // Examine the text in the response
                    response.text().then(function(data) {
                    if(extension !== '.html'){
                        data = '<pre><code>'+data+'</code></pre>';
                    }
                    else {
                        data = '<iframe srcdoc="'+sanitize_html(data)+'" style="border:none; width:100%; height:1000px;"></iframe>';
                    }
                    var btn_group = `<div class="btn-group" role="group">
                      <a class="btn btn-outline-secondary" href="${url}" target="_blank">Open in new tab</a>
                      <button type="button" class="btn btn-outline-secondary copy-url" data-target=${url}>Copy URL</button>
                      <button type="button" class="btn btn-outline-secondary copy-url" data-target=${s3_url}>Copy S3 URL</button>
                    </div>`;
                    var header = '<div class="card-header d-flex justify-content-between align-items-center">'+filename+btn_group+'</div>';
                    $("#file-preview").show();
                    $("#file-preview").html(header+'<div class="card-body">'+data+'</div>');
                });
            }
         )
         .catch(function(err) {
            console.log('Fetch Error :-S', err);
        });
    });

    //copy text to clipboard
    $('.toast').toast();

    $("body").on("click",".copy-url",function(e){
            var text = e.currentTarget.dataset.target;
            var $tmp = $("<input>");
            $("body").append($tmp);
            $tmp.val(text).select();
            document.execCommand("copy");
            $tmp.remove();
            $('#url-copied').toast('show');
        });


    // We are going to generate bucket/folder breadcrumbs. The resulting HTML will
    // look something like this:
    //
    // <li>Home</li>
    // <li>Library</li>
    // <li class="active">Samples</li>
    //
    // Note: this code is a little complex right now so it would be good to find
    // a simpler way to create the breadcrumbs.
    function folder2breadcrumbs(data) {
        // console.log('Bucket: ' + data.params.Bucket);
        // console.log('Prefix: ' + data.params.Prefix);

        if (data.params.Prefix && data.params.Prefix.length > 0) {
            // console.log('Set hash: ' + data.params.Prefix);
            window.location.hash = data.params.Prefix;
        } else {
            // console.log('Remove hash');
            removeHash();
        }

        // The parts array will contain the bucket name followed by all the
        // segments of the prefix, exploded out as separate strings.
        var parts = [data.params.Bucket];

        if (data.params.Prefix) {
            parts.push.apply(parts,
                data.params.Prefix.endsWith('/') ?
                data.params.Prefix.slice(0, -1).split('/') :
                data.params.Prefix.split('/'));
        }

        // console.log('Parts: ' + parts + ' (length=' + parts.length + ')');

        // Empty the current breadcrumb list
        $('#breadcrumb li').remove();

        // Now build the new breadcrumb list
        var buildprefix = '';
        $.each(parts, function(ii, part) {
            var ipart;

            // Add the bucket (the bucket is always first)
            if (ii === 0 || ii===1) {
                var a1 = $('<span>').text(part);
                ipart = $('<li>').addClass('breadcrumb-item').append(a1);
                // a1.click(function(e) {
                //     e.preventDefault();
                //     // console.log('Breadcrumb click bucket: ' + data.params.Bucket);
                //     s3exp_config = {
                //         Bucket: data.params.Bucket,
                //         Prefix: '',
                //         Delimiter: data.params.Delimiter
                //     };
                //     (s3exp_lister = s3list(s3exp_config, s3draw)).go();
                // });
                if(ii===1){
                    buildprefix += part + '/';
                }
                // Else add the folders within the bucket
            } else {
                buildprefix += part + '/';

                if (ii == parts.length - 1) {
                    ipart = $('<li>').addClass('breadcrumb-item active').text(part);
                } else {
                    var a2 = $('<a>').attr('href', '#').append(part);
                    ipart = $('<li>').addClass('breadcrumb-item').append(a2);

                    // Closure needed to enclose the saved S3 prefix
                    (function() {
                        var saveprefix = buildprefix;
                        // console.log('Part: ' + part + ' has buildprefix: ' + saveprefix);
                        a2.click(function(e) {
                            e.preventDefault();
                            console.log('Breadcrumb click object prefix: ' + saveprefix);
                            s3exp_config = {
                                Bucket: data.params.Bucket,
                                Prefix: saveprefix,
                                Delimiter: data.params.Delimiter
                            };
                            (s3exp_lister = s3list(s3exp_config, s3draw)).go();
                        });
                    })();
                }
            }
            $('#breadcrumb').append(ipart);
        });
    }

    function s3draw(data, complete) {
    $('li.li-bucket').remove();
    $("#file-preview").hide();
    folder2breadcrumbs(data);

    var path = data.params.Prefix.split("/")
    if(path.length > 3){
        $('#tb-s3objects').DataTable().rows.add([{
            Key: path.slice(0,path.length-2).join("/") + "/",
            render_name: false
        }]);
    }
    // Add each part of current path (S3 bucket plus folder hierarchy) into the breadcrumbs
    $.each(data.CommonPrefixes, function(i, prefix) {
        $('#tb-s3objects').DataTable().rows.add([{
            Key: prefix.Prefix,
            render_name: true
        }]);
    });

    // Add S3 objects to DataTable
    $('#tb-s3objects').DataTable().rows.add(data.Contents).draw();
    }

    function s3list(config, completecb) {
    // console.log('s3list config: ' + JSON.stringify(config));
    var params = {
        Bucket: config.Bucket,
        Prefix: config.Prefix,
        Delimiter: config.Delimiter
    };
    var scope = {
        Contents: [],
        CommonPrefixes: [],
        params: params,
        stop: false,
        completecb: completecb
    };

    return {
        // This is the callback that the S3 API makes when an S3 listObjectsV2
        // request completes (successfully or in error). Note that a single call
        // to listObjectsV2 may not be enough to get all objects so we need to
        // check if the returned data is truncated and, if so, make additional
        // requests with a 'next marker' until we have all the objects.
        cb: function(err, data) {
            if (err) {
                console.log('Error: ' + JSON.stringify(err));
                console.log('Error: ' + err.stack);
                scope.stop = true;
                $('#bucket-loader').removeClass('fa-spin');
                bootbox.alert("Error accessing S3 bucket " + scope.params.Bucket + ". Error: " + err);
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
                data.Contents = data.Contents.filter(function(el) {
                    return el.Key !== scope.params.Prefix;
                });

                // Optionally, filter the root index.html out of the listed S3 objects
                if (HIDE_INDEX) {
                    // console.log("Filter: remove index.html");
                    data.Contents = data.Contents.filter(function(el) {
                        return el.Key !== "index.html";
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
        go: function() {
            scope.cb = this.cb;
            $('#bucket-loader').addClass('fa-spin');
            $('#tb-s3objects').DataTable().clear();
            s3.makeUnauthenticatedRequest('listObjectsV2', scope.params, this.cb);
        },

        stop: function() {
            scope.stop = true;
            delete scope.params.ContinuationToken;
            if (scope.completecb) {
                scope.completecb(scope, false);
            }
            $('#bucket-loader').removeClass('fa-spin');
        }
    };
    }

    function promptForBucketInput() {
    bootbox.prompt("Please enter the S3 bucket name", function(result) {
        if (result !== null) {
            resetDepth();
            s3exp_config = {
                Bucket: result,
                Delimiter: '/'
            };
            (s3exp_lister = s3list(s3exp_config, s3draw)).go();
        }
    });
    }

    function resetDepth() {
    $('#tb-s3objects').DataTable().column(1).visible(false);
    $('input[name="optionsdepth"]').val(['folder']);
    $('input[name="optionsdepth"][value="bucket"]').parent().removeClass('active');
    $('input[name="optionsdepth"][value="folder"]').parent().addClass('active');
    }

    $(document).ready(function() {
    console.log('ready');

    // Click handler for refresh button (to invoke manual refresh)
    $('#bucket-loader').click(function(e) {
        if ($('#bucket-loader').hasClass('fa-spin')) {
            // To do: We need to stop the S3 list that's going on
            // bootbox.alert("Stop is not yet supported.");
            s3exp_lister.stop();
        } else {
            delete s3exp_config.ContinuationToken;
            (s3exp_lister = s3list(s3exp_config, s3draw)).go();
        }
    });

    // // Click handler for bucket button (to allow user to change bucket)
    // $('#bucket-chooser').click(function(e) {
    //     promptForBucketInput();
    // });

    // $('#hidefolders').click(function(e) {
    //     $('#tb-s3objects').DataTable().draw();
    // });

    // // Folder/Bucket radio button handler
    // $("input:radio[name='optionsdepth']").change(function() {
    //     // console.log("Folder/Bucket option change to " + $(this).val());
    //     // console.log("Change options: " + $("input[name='optionsdepth']:checked").val());
    //
    //     // If user selected deep then we do need to do a full list
    //     if ($(this).val() == 'bucket') {
    //         console.log("Switch to bucket");
    //         var choice = $(this).val();
    //         $('#tb-s3objects').DataTable().column(1).visible(choice === 'bucket');
    //         delete s3exp_config.ContinuationToken;
    //         delete s3exp_config.Prefix;
    //         s3exp_config.Delimiter = '';
    //         (s3exp_lister = s3list(s3exp_config, s3draw)).go();
    //         // Else user selected folder then can do a delimiter list
    //     } else {
    //         console.log("Switch to folder");
    //         $('#tb-s3objects').DataTable().column(1).visible(false);
    //         delete s3exp_config.ContinuationToken;
    //         delete s3exp_config.Prefix;
    //         s3exp_config.Delimiter = '/';
    //         (s3exp_lister = s3list(s3exp_config, s3draw)).go();
    //     }
    // });

    function renderObject(data, type, full) {
        if (isthisdocument(s3exp_config.Bucket, data)) {
            // console.log("is this document: " + data);
            return fullpath2filename(data);
        } else if (isfolder(data)) {
            // console.log("is folder: " + data);
            if(full.render_name){
                return '<i class="fad fa-folder"></i> <a data-s3="folder" data-prefix="' + sanitize_html(data) + '" href="' + object2hrefvirt(s3exp_config.Bucket, data) + '">' + prefix2folder(data) + '</a>';
            } else {
                return '<i class="fad fa-folder-open"></i> <a data-s3="folder" data-prefix="' + sanitize_html(data) + '" href="' + object2hrefvirt(s3exp_config.Bucket, data) + '">  ..</a>';
            }

        } else {
            var icon = '<i class="fad fa-file"></i> ';
            var preview_button = '';
            if(haspreview(full)){
                var preview_button = '<button class="btn btn-outline-secondary preview-file-btn ml-3"><i class="fad fa-file-search mr-1"></i> Preview file</button>';
            }
            // console.log("not folder/this document: " + data);
            return '<div class="d-flex justify-content-between align-items-center"><div>'+icon +'<a data-s3="object" href="' + object2hrefvirt(s3exp_config.Bucket, data) + '"download="' + fullpath2filename(data) + '">'+ fullpath2filename(data) + '</a></div>' + preview_button + '</div>';
        }
    }

    function renderFolder(data, type, full) {
        return isfolder(data) ? "" : fullpath2pathname(data);
    }

    // Initial DataTable settings
    $('#tb-s3objects').DataTable({
        info: false,
        paging: false,
        iDisplayLength: 50,
        order: [
            [1, 'asc'],
            [0, 'asc']
        ],
        aoColumnDefs: [{
            "aTargets": [0],
            "mData": "Key",
            "mRender": function(data, type, full) {
                return (type == 'display') ? renderObject(data, type, full) : data;
            },
            "sType": "key"
        }, {
            "aTargets": [1],
            "mData": "Key",
            "width": "50%",
            "mRender": function(data, type, full) {
                return renderFolder(data, type, full);
            }
        }, {
            "aTargets": [2],
            "mData": "LastModified",
            "width": "35%",
            "mRender": function(data, type, full) {
                return data ? moment(data).fromNow() : "";
            }
        }, {
            "aTargets": [3],
            "width": "15%",
            "mData": function(data, type, val) {
                return data.Size ? ((type == 'display') ? bytesToSize(data.Size) : data.Size) : "";
            }
        }, ]
    });
    $('#tb-s3objects').DataTable().column(s3exp_columns.key).visible(false);
    // console.log("jQuery version=" + $.fn.jquery);

    // Custom sort for the Key column so that folders appear before objects
    $.fn.dataTableExt.oSort['key-asc'] = function(a, b) {
        var x = (isfolder(a) ? "0-" + a : "1-" + a).toLowerCase();
        var y = (isfolder(b) ? "0-" + b : "1-" + b).toLowerCase();
        return ((x < y) ? -1 : ((x > y) ? 1 : 0));
    };

    $.fn.dataTableExt.oSort['key-desc'] = function(a, b) {
        var x = (isfolder(a) ? "1-" + a : "0-" + a).toLowerCase();
        var y = (isfolder(b) ? "1-" + b : "0-" + b).toLowerCase();
        return ((x < y) ? 1 : ((x > y) ? -1 : 0));
    };

    // // Allow user to hide folders
    // $.fn.dataTableExt.afnFiltering.push(function(oSettings, aData, iDataIndex) {
    //     // console.log("hide folders");
    //     return $('#hidefolders').is(':checked') ? !isfolder(aData[0]) : true;
    // });

    // Delegated event handler for S3 object/folder clicks. This is delegated
    // because the object/folder rows are added dynamically and we do not want
    // to have to assign click handlers to each and every row.
    $('#tb-s3objects').on('click', 'a', function(event) {
        event.preventDefault();
        var target = event.target;
        // console.log("target href=" + target.href);
        console.log("target dataset=" + JSON.stringify(target.dataset));

        // If the user has clicked on a folder then navigate into that folder
        if (target.dataset.s3 === "folder") {
            resetDepth();
            delete s3exp_config.ContinuationToken;
            s3exp_config.Prefix = target.dataset.prefix;
            // s3exp_config.Delimiter = $("input[name='optionsdepth']:checked").val() == "folder" ? "/" : "";
            s3exp_config.Delimiter = "/";
            (s3exp_lister = s3list(s3exp_config, s3draw)).go();
            // Else user has clicked on an object so download it in new window/tab
        } else {
            window.open(target.href, '_blank');
        }
        return false;
    });

    // Document URL typically looks like this for path-style URLs:
    // - https://s3.amazonaws.com/mybucket1/index.html
    // - https://s3-us-west-2.amazonaws.com/mybucket2/index.html
    //
    // Document URL typically looks like this for virtual-hosted-style URLs:
    // - https://mybucket1.s3.amazonaws.com/index.html
    // - https://mybucket2.s3-us-west-2.amazonaws.com/index.html
    //
    // Document URL typically looks like this for S3 website hosting:
    // - http://mybucket3.s3-website-us-east-1.amazonaws.com/
    // - http://mybucket4.s3-website.eu-central-1.amazonaws.com/

    // TODO: need to support S3 website hosting option
    //
    // If we're launched from a bucket then let's try to determine the bucket
    // name so we can query it immediately, without requiring the user to
    // supply the bucket name.
    //
    // If the region was anything other than US Standard then we will also need
    // to infer the region so that we can initialize the S3 SDK properly.
    // console.log("Document URL: " + document.URL);
    // var urls = document.URL.split('/');
    // console.log("URL split: " + urls);

    // Using technique from https://gist.github.com/jlong/2428561
    // to parse the document URL.
    // var parser = document.createElement('a');
    // parser.href = document.URL;

    // URL format is scheme://[user:password@]domain:port/path?query_string#fragment_id
    // For example: http://example.com:3000/path/?name=abc#topic
    // console.log("protocol: " + parser.protocol); // => "http:"
    // console.log("hostname: " + parser.hostname); // => "example.com"
    // console.log("port    : " + parser.port); // => "3000"
    // console.log("pathname: " + parser.pathname); // => "/path/"
    // console.log("search  : " + parser.search); // => "?name=abc"
    // console.log("hash    : " + parser.hash); // => "#topic"
    // console.log("host    : " + parser.host); // => "example.com:3000"

    // If initial bucket has been hard-coded above then use it, else try to
    // derive the initial bucket from the document URL (useful if index.html was
    // launched directly from within a bucket), else prompt the user.
    if (s3exp_config.Bucket) {
        (s3exp_lister = s3list(s3exp_config, s3draw)).go(); //initial load
            }
    // else if (parser.hostname.endsWith('amazonaws.com')) {
    //             // Hostname is likely to be in one of the following forms:
    //             // - s3.amazonaws.com
    //             // - bucket1.s3.amazonaws.com
    //             // - s3-us-west-2.amazonaws.com
    //             // - bucket2.s3-us-west-2.amazonaws.com
    //
    //             // If using static website hosting, the hostname will be of the form:
    //             // - bucket.s3-website.eu-central-1.amazonaws.com/
    //
    //             // The following is also a legal form, but we do not support it, so
    //             // you should use s3-eu-central-1 rather than s3.eu-central-1:
    //             // - bucket3.s3.eu-central-1.amazonaws.com
    //
    //             var bucket;
    //             var region;
    //             var hostnames = parser.hostname.split('.');
    //             var pathnames = parser.pathname.split('/');
    //
    // //             console.log("count of words in hostname=" + hostnames.length);
    // //             console.log("count of words in pathname=" + pathnames.length);
    // //
    // //             console.log("hostnames=" + hostnames);
    // //             console.log("pathnames=" + pathnames);
    //
    //             // If bucket prefix not included in hostname
    //             if (hostnames[0].match(/^s3-/) || hostnames[0].match(/^s3$/)) {
    //                 bucket = pathnames[1];
    //                 region = hostnames[0];
    //                 // console.log("path bucket=" + bucket);
    //                 // console.log("path region=" + region);
    //             } else {
    //                 bucket = hostnames[0];
    //                 region = hostnames[hostnames.length - 3];
    //                 // console.log("host bucket=" + bucket);
    //                 // console.log("host region=" + region);
    //             }
    //
    //             // If we found statically-hosted website prefix or explicit region, for
    //             // example s3-us-west-2, then get region else use the default of US Standard
    //             if (region !== 's3') {
    //                 if (region.startsWith('s3-website-')) {
    //                     AWS.config.region = region.substring(11);
    //                 } else if (region.startsWith('s3-') || region.startsWith('s3.')) {
    //                     AWS.config.region = region.substring(3);
    //                 } else {
    //                     AWS.config.region = region;
    //                 }
    //             }
    //
    //             // console.log("AWS region=" + AWS.config.region);
    //             // console.log("S3 bucket=" + bucket);
    //
    //             // Create and initialize S3 object
    //             s3 = new AWS.S3();
    //             s3exp_config = {
    //                 Bucket: bucket,
    //                 Delimiter: '/'
    //             };
    //
    //             if (window.location.hash) {
    //                 // console.log("Location hash=" + window.location.hash);
    //                 s3exp_config.Prefix = window.location.hash.substring(1);
    //             }
    //
    //             // Do initial bucket list
                // (s3exp_lister = s3list(s3exp_config, s3draw)).go();
    //         } else {
    //             promptForBucketInput();
    //         }
    });
});
