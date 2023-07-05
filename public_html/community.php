<?php
// ini_set('display_errors', 1);
// ini_set('display_startup_errors', 1);
// error_reporting(E_ALL);

$title = 'Community';
$subtitle = 'Find out who is involved in the sanger-tol project';
$md_github_url = 'https://github.com/sanger-tol/pipelines-website/blob/updating_stats/sanger-tol-partners.yaml';
$import_leaflet = true;
include '../includes/header.php';
?>

<h1>Introduction</h1>
<p>sanger-tol is by design a collaborative effort, and would not exist if it were not for the efforts of many dedicated contributors.</p>
<ul>
    <li><a href="#contributors">Contributors</a></li>
    <li><a href="#organisations">Organisations</a></li>
    <li><a href="#initiatives">Projects we are involved with</a></li>
</ul>

<?php echo _h1('Contributors'); ?>
<p>The sanger-tol pipelines and community is driven by many individuals, listed below. This list updates automatically.</p>
<p>Want to see who's working on what? See the <a href="/stats#contributor_leaderboard">contributor leaderboard</a> on the Statistics page.</p>
<p class="pt-3">
    <?php
    $stats_json_fn = dirname(dirname(__FILE__)) . '/nfcore_stats.json';
    $stats_json = json_decode(file_get_contents($stats_json_fn));
    $contributors = [];
    foreach (['pipelines', 'core_repos'] as $repo_type) {
        foreach ($stats_json->{$repo_type} as $repo) {
            foreach ($repo->contributors as $contributor) {
                $contributors[$contributor->author->login] = $contributor->author;
            }
        }
    }
    // Random order!
    $logins = array_keys($contributors);
    shuffle($logins);
    foreach ($logins as $login) {
        $author = $contributors[$login];
        echo '<a title="@' .
            $author->login .
            '" href="' .
            $author->html_url .
            '" target="_blank" data-bs-toggle="tooltip"><img src="' .
            $author->avatar_url .
            '" class="border rounded-circle me-1 mb-1" width="50" height="50"></a>';
    }
    ?>
</p>

<?php echo _h1('Partners'); ?>
<p>Tree of Life programme at Wellcome Sanger Institute works with several partners across the globe to deliver its goal to sequence all eukaryotic life on the planet. Some of these organisations are listed below, along with a key person who you can contact for advice.</p>
<blockquote>Is your group missing? Please submit a pull request to add yourself! It's just a few lines in a <a href="https://github.com/sanger-tol/pipelines-website/blob/main/sanger-tol-partners.yaml">simple YAML file.</a></blockquote>

<div class="card contributors-map-card">
    <div class="card-body" id="contributors-map"></div>
</div>
<div class="row row-cols-1 row-cols-md-2 g-4">

    <?php
    // Parse YAML contributors file
    require '../vendor/autoload.php';

    use Spyc;

    $locations = [];
    $contributors = spyc_load_file('../sanger-tol-partners.yaml');
    $contributors_html = '';
    foreach ($contributors['contributors'] as $idx => $c) {
        // Start card div
        $contributors_html .= '<div class="col"><div class="card contributor h-100"><div class="card-body">';
        // Header, title
        $img_path = '';
        if (array_key_exists('image_fn', $c)) {
            // Dark theme
            $hide_dark = '';
            $dark_img_path = 'assets/img/contributors-white/' . $c['image_fn'];
            if ($c['image_fn'] and file_exists($dark_img_path)) {
                $contributors_html .=
                    '<img class="contributor_logo hide-light" title="' .
                    $c['full_name'] .
                    '" src="' .
                    $dark_img_path .
                    '">';
                $hide_dark = 'hide-dark';
            }
            // Normal, light theme
            $img_path = 'assets/img/contributors-colour/' . $c['image_fn'];
            if ($c['image_fn'] and file_exists($img_path)) {
                $contributors_html .=
                    '<img class="contributor_logo ' .
                    $hide_dark .
                    '" title="' .
                    $c['full_name'] .
                    '" src="' .
                    $img_path .
                    '">';
            } else {
                $img_path = '';
            }
        }
        $card_id = $idx;
        if (array_key_exists('full_name', $c)) {
            $card_id = preg_replace('/[^a-z]+/', '-', strtolower($c['full_name']));
            $contributors_html .= '<h5 class="card-title" id="' . $card_id . '">';
            if (array_key_exists('url', $c)) {
                $contributors_html .= ' <a href="' . $c['url'] . '" target="_blank">';
            }
            $contributors_html .= $c['full_name'];
            if (array_key_exists('url', $c)) {
                $contributors_html .= ' </a>';
            }
            $contributors_html .= '</h5>';
        }
        if (array_key_exists('affiliation', $c)) {
            $contributors_html .= '<h6 class="card-subtitle mb-2 text-muted">';
            if (array_key_exists('affiliation_url', $c)) {
                $contributors_html .= '<a href="' . $c['affiliation_url'] . '" target="_blank">';
            }
            $contributors_html .= $c['affiliation'];
            if (array_key_exists('affiliation_url', $c)) {
                $contributors_html .= '</a>';
            }
            $contributors_html .= '</h6>';
        }
        // Description
        if (array_key_exists('description', $c)) {
            $contributors_html .= '<p class="card-text small text-muted">' . $c['description'] . '</p> ';
        }
        // Contact person
        $contributors_html .= '<div class="contributor_contact">';
        if (array_key_exists('contact_email', $c)) {
            $contributors_html .=
                '<a href="mailto:' .
                $c['contact_email'] .
                '" class="badge bg-light text-dark fw-normal" data-bs-toggle="tooltip" title="Primary contact: ' .
                $c['contact_email'] .
                '"><i class="far fa-envelope"></i> ';
            if (array_key_exists('contact', $c)) {
                $contributors_html .= $c['contact'];
            } else {
                $contributors_html .= $c['contact_email'];
            }
            $contributors_html .= '</a> ';
        } elseif (array_key_exists('contact', $c)) {
            $contributors_html .= '<span class="badge bg-light text-dark fw-normal">' . $c['contact'] . '</span> ';
        }
        if (array_key_exists('contact_github', $c)) {
            $contributors_html .=
                '<a href="https://github.com/' .
                trim($c['contact_github'], '@') .
                '/" target="_blank" class="badge bg-light text-dark fw-normal" data-bs-toggle="tooltip" title="Primary contact: GitHub @' .
                trim($c['contact_github'], '@') .
                '"><i class="fab fa-github"></i> ' .
                trim($c['contact_github'], '@') .
                '</a> ';
        }
        if (array_key_exists('twitter', $c)) {
            $contributors_html .=
                '<a href="https://twitter.com/' .
                trim($c['twitter'], '@') .
                '/" target="_blank" class="badge bg-light text-dark fw-normal" data-bs-toggle="tooltip" title="Institutional twitter: @' .
                trim($c['twitter'], '@') .
                '"><i class="fab fa-twitter"></i> @' .
                trim($c['twitter'], '@') .
                '</a> ';
        }
        $contributors_html .= '</div>';
        // Close card div
        $contributors_html .= '</div></div></div>';

        // Location JSON
        if (array_key_exists('location', $c)) {
            $location['location'] = $c['location'];
            $location['full_name'] = array_key_exists('full_name', $c) ? $c['full_name'] : '';
            $location['card_id'] = $card_id;
            if ($img_path) {
                $location['image'] =
                    '<br><a href="#' .
                    $card_id .
                    '"><img class="contributor_map_logo" title="' .
                    $location['full_name'] .
                    '" src="' .
                    $img_path .
                    '"></a>';
            } else {
                $location['image'] = '';
            }
            array_push($locations, $location);
        }
    }

    echo $contributors_html;
    ?>

</div>
<script type="text/javascript">
    var locations = <?php echo json_encode($locations, JSON_PRETTY_PRINT); ?>;

    $(function() {
        var map = L.map('contributors-map', {
            zoom: 2
        });
        var greenIcon = new L.Icon({
            iconUrl: 'assets/img/marker-icon-2x-green.png',
            shadowUrl: 'assets/img/marker-shadow.png',
            iconSize: [25, 41],
            iconAnchor: [12, 41],
            popupAnchor: [1, -34],
            shadowSize: [41, 41]
        });

        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }).addTo(map);

        var latlngs = [];
        locations.forEach(function(marker) {
            if (marker != null) {
                L.marker(marker.location, {
                    icon: greenIcon
                }).addTo(map).bindPopup('<a href="#' + marker.card_id + '">' + marker.full_name + '</a>' + marker.image);
                latlngs.push(marker.location);
            }
        });
        map.fitBounds(latlngs);
    });
</script>

<h1 id="initiatives">Projects we are involved with<a href="#initiatives" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a></h1>
<p>A number of large projects are committing to working with the nf-core community for data analysis pipelines. You can see an overview of these below.</p>

<h3 id="dfg_testimonial">
    <img width="350px" src="/assets/img/contributors-colour/dfg_logo.svg" class="float-end ps-4" />

    German National Sequencing Initiative
    <a href="#dfg_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
</h3>
<p><a href="https://www.dfg.de/en/service/press/press_releases/2018/press_release_no_06/index.html" target="_blank">The German Funding Body (DFG)</a>
    has approved funding to establish 4 national high-throughput sequencing centers in Germany. The project will rely on <em>nf-core</em> pipelines for analyzing
    large-scale genomics data. Contributors from the Kiel, Tübingen and Dresden sites are already actively contributing to nf-core, and the other sequencing center
    in Cologne/Bonn (West German Genome Center) has committed to joining and contributing its expertise.</p>
<div class="clearfix"></div>

<h3 id="easi_genomics_testimonial">
    <img width="350px" src="/assets/img/contributors-<?php echo $theme == 'dark'
        ? 'white'
        : 'colour'; ?>/EASI-Genomics.svg" class="float-end ps-4 darkmode-image" />
    EASI-Genomics
    <a href="#easi_genomics_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
</h3>
<p><a href="https://www.easi-genomics.eu/" target="_blank">EASI-Genomics</a>
    provides easy and seamless access to cutting-edge DNA sequencing technologies to researchers from academia and industry, within a framework that ensures compliance with ethical and legal requirements, as well as FAIR and secure data management.
    EASI-Genomics aims to build a community of practice, which leverages advanced sequencing technologies beyond country and sector borders to tackle global challenges in science.
    It organizes competitive calls to get fully-funded access to state-of-the-art facilities and services, including integrative data analysis.</p>
<p>EASI-Genomics partners are committed to working within community best-practices for the bioinformatics processing of the data that they produce.
    To this end they will use <em>nf-core</em> pipelines to process data and contribute new developments and pipelines within the nf-core framework.</p>
<div class="clearfix"></div>

<h3 id="bovreg_testimonial">
    <img width="350px" src="/assets/img/contributors-colour/bovreg.svg" class="float-end ps-4" />
    BovReg
    <a href="#bovreg_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
</h3>
<p><a href="https://www.bovreg.eu/" target="_blank">BovReg</a> is a
    <a href="https://cordis.europa.eu/project/id/815668" target="_blank">H2020 project</a>
    to form a consortium that will provide a comprehensive map of functionally active genomic features in cattle and how their (epi)genetic variation in beef and dairy breeds translates into phenotypes.
</p>
<p>BovReg reference bioinformatics pipelines will adhere to <em>nf-core</em> guidelines in order to deliver reproducible data analyses to the community.</p>
<div class="clearfix"></div>
<h3 id="dockstore">
    <img width="350px" src="/assets/img/contributors-colour/dockstore.svg" class="float-end ps-0" />
    Dockstore
    <a href="#dockstore" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
</h3>
<p><a href="https://www.dockstore.org/" target="_blank">Dockstore</a> is a platform used by researchers across the world to share reproducible computational analyses and hosts hundreds of individual tools and workflows created by more than a hundred different contributors.</p>
<p>By using lightweight containerization technology along with the essential metadata needed for combining tools into scientific analysis "recipes", Dockstore allows users to create, share, publish (using citable DOIs) and reproducibly reuse these recipes across platforms and compute environments.</p>
<div class="clearfix"></div>

<h3 id="workflowhub">
    <img width="250px" src="/assets/img/contributors-colour/workflowhub.svg" class="float-end px-4" />
    Workflow Hub
    <a href="#workflowhub" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
</h3>
<p><a href="https://www.workflowhub.eu/" target="_blank">Workflow Hub</a> was created as part of the <a href="https://www.eosc-life.eu/" target="_blank">EOSC-Life</a> WP2: <em>Tools Collaboratory and Research Objects</em> to glue in federated workflow and tool descriptions across research infrastructures.</p>
<p>It is workflow system agnostic, supports a repository of workflows in native and standardised form and the virtual aggregation of established tool, workflow and registries to support discovery over a fragmented ecosystem. The federated registry will support a common API to simplify access for tool developers.</p>
<div class="clearfix"></div>

<?php include '../includes/footer.php';
