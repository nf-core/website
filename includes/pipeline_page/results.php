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
                    <!-- Title and breadcrumbs, wrap them on small displays and reshuffle order to put breadcrumb last -->
                    <div class="d-flex align-items-center title-bar flex-wrap flex-md-nowrap">
                        <!-- App title -->
                        <div class="title order-0 order-md-0">
                            <i class="fa-regular fa-folder-tree fa-lg me-md-3"></i>
                        </div>
                        <!-- Bucket breadcrumbs -->
                        <ul id="breadcrumb" class="breadcrumb order-1 order-md-0 mb-0">
                            <li class="breadcrumb-item active">
                                <a href="#"><i class="fas fa-circle-notch fa-spin"></i></a>
                            </li>
                        </ul>
                        <div class="ms-auto order-0 order-md-0">
                            <button type="button" class="btn btn-outline-secondary copy-url ms-3" data-bs-target="">Copy Bucket S3 URL</button>
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

<div class="toast" id="url-copied" role="alert" aria-live="assertive" aria-atomic="true">
    <div class="toast-header">
        <img src="/assets/img/logo/nf-core-logo-square.png" class="rounded me-2" alt="">
        <strong class="me-auto">URL copied to clipboard!</strong>
        <button type="button" class="ms-2 mb-1 btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
</div>

<script type="text/javascript">
    var HIDE_INDEX = true;
    var prefix = '<?php echo $pipeline->name; ?>/results-<?php echo $release_hash; ?>/';
    var suffix = '';
    if (window.location.hash.length > 0) {
        if (window.location.hash.substr(1).split("/")[0] === "<?php echo $pipeline->name; ?>") {
            prefix = window.location.hash.substr(1);
            if (!prefix.endsWith('/')) {
                suffix = prefix.split("/").pop();
                prefix = prefix.substr(0, prefix.lastIndexOf("/") + 1);
            }

        }
    }
    var s3exp_config = {
        Region: 'us-east-1',
        Bucket: 'sanger-tol-tests',
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
