<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.
$mainpage_container = false;
$aws = true;
ob_start();
?>
<div id="page-wrapper">
    <div class="row">
        <div class="col-lg-12">
            <div class="row ">
                <!-- Panel including title, breadcrumbs, and controls -->
                <div class="col">
                    <!-- Title and breadcrumbs -->
                    <div class="d-flex align-items-center title-bar">
                        <!-- App title -->
                        <div class="title ">
                            <i class="fab fa-aws fa-lg mr-3 mb-3"></i>
                        </div>
                        <!-- Bucket breadcrumbs -->
                        <ul id="breadcrumb" class="breadcrumb">
                            <li class="breadcrumb-item active">
                                <a href="#"><i class="fas fa-circle-notch fa-spin"></i></a>
                            </li>
                        </ul>
                        <button type="button" class="btn btn-outline-secondary copy-url ml-auto" data-target="">Copy Bucket S3 URL</button>
                    </div>
                    <!-- Controls -->
                    <div id="navbuttons" class="">
                        <div>
                            <!-- Dual purpose: progress spinner and refresh button, plus object count -->
                            <div class="btn-group" id="refresh">
                                <span id="bucket-loader" class="btn fa fa-refresh fa-2x " title="Refresh"></span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <!-- Panel including S3 object table -->
            <div class="panel-body">
                <table class="table table-bordered table-hover responsive" id="tb-s3objects">
                    <thead>
                        <tr>
                            <th>Name</th>
                            <th>Folder</th>
                            <th class="text-left">Last Modified</th>
                            <th class="text-left">Size</th>
                        </tr>
                    </thead>
                    <tbody id="tbody-s3objects"></tbody>
                </table>
            </div>
        </div>
    </div>
</div>
<div id="file-preview" class="card">
</div>

<div class="toast" id="url-copied" data-delay="5000" role="alert" aria-live="assertive" aria-atomic="true">
    <div class="toast-header">
        <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded mr-2" alt="">
        <strong class="mr-auto">URL copied to clipboard!</strong>
        <button type="button" class="ml-2 mb-1 close" data-dismiss="toast" aria-label="Close">
            <span aria-hidden="true">&times;</span>
        </button>
    </div>
</div>

<script type="text/javascript">
    var HIDE_INDEX = true;
    var prefix = '<?php echo $pipeline->name ?>/results-<?php echo $release_hash ?>/';
    var suffix = '';
    if (window.location.hash.length > 0) {
        if (window.location.hash.substr(1).split("/")[0] === "<?php echo $pipeline->name ?>") {
            prefix = window.location.hash.substr(1);
            if (!prefix.endsWith('/')) {
                suffix = prefix.split("/").pop();
                prefix = prefix.substr(0, prefix.lastIndexOf("/") + 1);
            }

        }
    }
    var s3exp_config = {
        Region: 'eu-west-1',
        Bucket: 'nf-core-awsmegatests',
        Prefix: prefix,
        Suffix: suffix,
        Delimiter: '/'
    };
    var s3exp_lister = null;
    var s3exp_columns = {
        key: 1,
        folder: 2,
        date: 3,
        size: 4
    };
</script>

<?php
$content = ob_get_contents();
ob_end_clean();
?>