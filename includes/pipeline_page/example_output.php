<?php
// Build the HTML for a pipeline documentation page.
// Imported by public_html/pipeline.php - pulls a markdown file from GitHub and renders.
$mainpage_container = false;
$aws=true;
ob_start();
?>
<div class="container container-xl">
    <div id="page-wrapper">
        <div class="row">
            <div class="col-lg-12">
                <div class="row ">
                    <!-- Panel including title, breadcrumbs, and controls -->
                    <div class="col">
                        <!-- Title and breadcrumbs -->
                        <div class="d-flex align-items-center">
                            <!-- App title -->
                            <div class="title ">
                                <i class="fab fa-aws fa-lg mr-3 mb-3"></i>
                            </div>
                            <!-- Bucket breadcrumbs -->
                            <div class="">
                                <ul id="breadcrumb" class="breadcrumb">
                                    <li class="breadcrumb-item active">
                                        <a href="#"><i class="fas fa-circle-notch fa-spin"></i></a>
                                    </li>
                                </ul>
                            </div>
                        </div>
                        <!-- Controls -->
                        <div id="navbuttons" class="">
                            <div>
                                <!-- Dual purpose: progress spinner and refresh button, plus object count -->
                                <div class="btn-group" id="refresh">
                                    <span id="bucket-loader" class="btn fa fa-refresh fa-2x " title="Refresh"></span>
                                    <!-- <span id="badgecount" class="badge pull-right" title="Object count">42</span> -->
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
                                <th>Last Modified</th>
                                <th>Size</th>
                            </tr>
                        </thead>
                        <tbody id="tbody-s3objects"></tbody>
                    </table>
                </div>
                </div>
            </div>
        </div>
    </div>
    <div id="file-preview" class="card">
    </div>
</div>
<script type="text/javascript">
    var HIDE_INDEX = true;
    var s3exp_config = {
        Region: 'eu-west-1',
        Bucket: 'nf-core-awsmegatests',
        Prefix: '<?php echo $pipeline->name?>/results-<?php echo $release_hash ?>/',
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
