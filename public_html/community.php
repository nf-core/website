<?php
// ini_set('display_errors', 1);
// ini_set('display_startup_errors', 1);
// error_reporting(E_ALL);

$title = 'Community';
$subtitle = 'Find out who is involved in the nf-core project';
$md_github_url = 'https://github.com/nf-core/nf-co.re/blob/master/nf-core-contributors.yaml';
$import_leaflet = true;
include('../includes/header.php');

?>

<h1>Introduction</h1>
<p>nf-core is by design a collaborative effort, and would not exist if it were not for the efforts of many dedicated contributors.</p>
<ul>
    <li><a href="#contributors">Contributors</a></li>
    <li><a href="#organisations">Organisations</a></li>
    <li><a href="#initiatives">Projects we are involved with</a></li>
    <li><a href="#testimonials">Testimonials</a></li>
</ul>


<h1 id="contributors"><a href="#contributors" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Contributors</h1>
<p>The nf-core pipelines and community is driven by many individuals, listed below. This list updates automatically.</p>
<p>Want to see who's working on what? See the <a href="/stats#contributor_leaderboard">contributor leaderboard</a> on the Statistics page.</p>
<p class="pt-3">
<?php
$stats_json_fn = dirname(dirname(__FILE__)).'/nfcore_stats.json';
$stats_json = json_decode(file_get_contents($stats_json_fn));
$contributors = [];
foreach(['pipelines', 'core_repos'] as $repo_type){
    foreach($stats_json->{$repo_type} as $repo){
        foreach($repo->contributors as $contributor){
            $contributors[$contributor->author->login] = $contributor->author;
        }
    }
}
// Random order!
$logins = array_keys($contributors);
shuffle($logins);
foreach($logins as $login){
    $author = $contributors[$login];
    echo '<a title="@'.$author->login.'" href="'.$author->html_url.'" target="_blank" data-toggle="tooltip"><img src="'.$author->avatar_url.'" class="border rounded-circle mr-1 mb-1" width="50" height="50"></a>';
}
?>
</p>

<h1 id="organisations"><a href="#organisations" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Organisations</h1>
<p>Some of the organisations running nf-core pipelines are listed below, along with a key person who you can contact for advice.</p>
<blockquote>Is your group missing? Please submit a pull request to add yourself! It's just a few lines in a <a href="https://github.com/nf-core/nf-co.re/blob/master/nf-core-contributors.yaml">simple YAML file.</a></blockquote>

<div class="card contributors-map-card">
    <div class="card-body" id="contributors-map"></div>
</div>
<div class="card-deck">

<?php
// Parse YAML contributors file
$locations = [];
require_once("../includes/libraries/Spyc.php");
$contributors = spyc_load_file('../nf-core-contributors.yaml');
$contributors_html = '';
foreach($contributors['contributors'] as $idx => $c){
    // Start card div
    $contributors_html .= '<div class="card contributor card_deck_card"><div class="card-body">';
    // Header, title
    $img_path = '';
    if(array_key_exists('image_fn', $c)){
        // Dark theme
        $hide_dark = '';
        $dark_img_path = 'assets/img/contributors-white/'.$c['image_fn'];
        if($c['image_fn'] and file_exists($dark_img_path)){
            $contributors_html .= '<img class="contributor_logo hide-light" title="'.$c['full_name'].'" src="'.$dark_img_path.'">';
            $hide_dark = 'hide-dark';
        }
        // Normal, light theme
        $img_path = 'assets/img/contributors-colour/'.$c['image_fn'];
        if($c['image_fn'] and file_exists($img_path))
            $contributors_html .= '<img class="contributor_logo '.$hide_dark.'" title="'.$c['full_name'].'" src="'.$img_path.'">';
        else $img_path = '';
    }
    $card_id = $idx;
    if(array_key_exists('full_name', $c)){
        $card_id = preg_replace('/[^a-z]+/', '-', strtolower($c['full_name']));
        $contributors_html .= '<h5 class="card-title" id="'.$card_id.'">';
        if(array_key_exists('url', $c))
            $contributors_html .= ' <a href="'.$c['url'].'" target="_blank">';
        $contributors_html .= $c['full_name'];
        if(array_key_exists('url', $c))
            $contributors_html .= ' </a>';
        $contributors_html .= '</h5>';
    }
    if(array_key_exists('affiliation', $c)){
        $contributors_html .= '<h6 class="card-subtitle mb-2 text-muted">';
        if(array_key_exists('affiliation_url', $c))
            $contributors_html .= '<a href="'.$c['affiliation_url'].'" target="_blank">';
        $contributors_html .= $c['affiliation'];
        if(array_key_exists('affiliation_url', $c))
            $contributors_html .= '</a>';
        $contributors_html .= '</h6>';
    }
    // Description
    if(array_key_exists('description', $c))
        $contributors_html .= '<p class="small text-muted">'.$c['description'].'</p> ';
    // Contact person
    $contributors_html .= '<div class="contributor_contact">';
    if(array_key_exists('contact_email', $c)){
        $contributors_html .= '<a href="mailto:'.$c['contact_email'].'" class="badge badge-light" data-toggle="tooltip" title="Primary contact: '.$c['contact_email'].'"><i class="far fa-envelope"></i> ';
        if(array_key_exists('contact', $c))
            $contributors_html .= $c['contact'];
        else
            $contributors_html .= $c['contact_email'];
        $contributors_html .= '</a> ';
    }
    else if(array_key_exists('contact', $c))
        $contributors_html .= '<span class="badge badge-light">'.$c['contact'].'</span> ';
    if(array_key_exists('contact_github', $c))
        $contributors_html .= '<a href="https://github.com/'.trim($c['contact_github'], '@').'/" target="_blank" class="badge badge-light" data-toggle="tooltip" title="Primary contact: GitHub @'.trim($c['contact_github'], '@').'"><i class="fab fa-github"></i> '.trim($c['contact_github'], '@').'</a> ';
    if(array_key_exists('twitter', $c))
        $contributors_html .= '<a href="https://twitter.com/'.trim($c['twitter'], '@').'/" target="_blank" class="badge badge-light" data-toggle="tooltip" title="Institutional twitter: @'.trim($c['twitter'], '@').'"><i class="fab fa-twitter"></i> @'.trim($c['twitter'], '@').'</a> ';
    $contributors_html .= '</div>';
    // Close card div
    $contributors_html .= '</div></div>';

    // Location JSON
    if(array_key_exists('location', $c)){
        $location['location'] = $c['location'];
        $location['full_name'] = array_key_exists('full_name', $c) ? $c['full_name'] : '';
        $location['card_id'] = $card_id;
        if($img_path){
            $location['image'] = '<br><a href="#'.$card_id.'"><img class="contributor_map_logo" title="'.$location['full_name'].'" src="'.$img_path.'"></a>';
        } else $location['image'] = '';
        array_push($locations, $location);
    }
}

echo $contributors_html;
?>

</div>
<script type="text/javascript">
var locations = <?php echo json_encode($locations, JSON_PRETTY_PRINT); ?>;

$(function(){
    var map = L.map('contributors-map', {
        fullscreenControl: true,
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
            L.marker(marker.location, {icon: greenIcon}).addTo(map).bindPopup('<a href="#'+marker.card_id+'">'+marker.full_name+'</a>'+marker.image);
            latlngs.push(marker.location);
        }
    });
    map.fitBounds(latlngs);
});
</script>

<h1 id="initiatives"><a href="#initiatives" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Projects we are involved with</h1>
<p>A number of large projects are committing to working with the nf-core community for data analysis pipelines. You can see an overview of these below.</p>

<h3 id="dfg_testimonial"><a href="#dfg_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
    <img width="350px" src="/assets/img/contributors-colour/dfg_logo.svg" class="float-right pl-4" />
    German National Sequencing Initiative
</h3>
<p><a href="https://www.dfg.de/en/service/press/press_releases/2018/press_release_no_06/index.html" target="_blank">The German Funding Body (DFG)</a>
has approved funding to establish 4 national high-throughput sequencing centers in Germany. The project will rely on <em>nf-core</em> pipelines for analyzing
large-scale genomics data. Contributors from the Kiel, Tübingen and Dresden sites are already actively contributing to nf-core, and the other sequencing center
in Cologne/Bonn (West German Genome Center) has committed to joining and contributing its expertise.</p>
<div class="clearfix"></div>

<h3 id="easi_genomics_testimonial"><a href="#easi_genomics_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
    <img width="350px" src="/assets/img/contributors-<?php echo $theme == 'dark' ? 'white' : 'colour'; ?>/EASI-Genomics.svg" class="float-right pl-4 darkmode-image" />
    EASI-Genomics
</h3>
<p><a href="https://www.easi-genomics.eu/" target="_blank">EASI-Genomics</a>
provides easy and seamless access to cutting-edge DNA sequencing technologies to researchers from academia and industry, within a framework that ensures compliance with ethical and legal requirements, as well as FAIR and secure data management.
EASI-Genomics aims to build a community of practice, which leverages advanced sequencing technologies beyond country and sector borders to tackle global challenges in science.
It organizes competitive calls to get fully-funded access to state-of-the-art facilities and services, including integrative data analysis.</p>
<p>EASI-Genomics partners are committed to working within community best-practices for the bioinformatics processing of the data that they produce.
To this end they will use <em>nf-core</em> pipelines to process data and contribute new developments and pipelines within the nf-core framework.</p>
<div class="clearfix"></div>

<h3 id="bovreg_testimonial"><a href="#bovreg_testimonial" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
    <img width="350px" src="/assets/img/contributors-colour/bovreg.svg" class="float-right pl-4" />
    BovReg
</h3>
<p><a href="https://www.bovreg.eu/" target="_blank">BovReg</a> is a
<a href="https://cordis.europa.eu/project/id/815668" target="_blank">H2020 project</a>
to form a consortium that will provide a comprehensive map of functionally active genomic features in cattle and how their (epi)genetic variation in beef and dairy breeds translates into phenotypes.</p>
<p>BovReg reference bioinformatics pipelines will adhere to <em>nf-core</em> guidelines in order to deliver reproducible data analyses to the community.</p>
<div class="clearfix"></div>
<h3 id="dockstore"><a href="#dockstore" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
    <img width="350px" src="/assets/img/contributors-colour/dockstore.svg" class="float-right pl-0" />
    Dockstore
</h3>
<p><a href="https://www.dockstore.org/" target="_blank">Dockstore</a> is a platform used by researchers across the world to share reproducible computational analyses and hosts hundreds of individual tools and workflows created by more than a hundred different contributors.</p>
<p>By using lightweight containerization technology along with the essential metadata needed for combining tools into scientific analysis "recipes", Dockstore allows users to create, share, publish (using citable DOIs) and reproducibly reuse these recipes across platforms and compute environments.</p>
<div class="clearfix"></div>

<h3 id="workflowhub"><a href="#workflowhub" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>
    <img width="250px" src="/assets/img/contributors-colour/workflowhub.svg" class="float-right px-4" />
    Workflow Hub
</h3>
<p><a href="https://www.workflowhub.eu/" target="_blank">Workflow Hub</a> was created as part of the <a href="https://www.eosc-life.eu/" target="_blank">EOSC-Life</a> WP2: <em>Tools Collaboratory and Research Objects</em> to glue in federated workflow and tool descriptions across research infrastructures.</p>
<p>It is workflow system agnostic, supports a repository of workflows in native and standardised form and the virtual aggregation of established tool, workflow and registries to support discovery over a fragmented ecosystem. The federated registry will support a common API to simplify access for tool developers.</p>
<div class="clearfix"></div>

<h1 id="testimonials"><a href="#testimonials" class="header-link"><span class="fas fa-link" aria-hidden="true"></span></a>Testimonials</h1>

<div class="card-columns">

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">.<a href="https://twitter.com/nextflowio?ref_src=twsrc%5Etfw">@nextflowio</a> and <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a> are making me hate bioinformatics a bit less</p>&mdash; Hadrien Gourlé (@HGourle) <a href="https://twitter.com/HGourle/status/1158973454160465920?ref_src=twsrc%5Etfw">August 7, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">Just tested <a href="https://t.co/wERtIS7xPp">https://t.co/wERtIS7xPp</a> with bacterial WGS data. It was very easy to work with. I am very impressed by nf-core! Great community and very important project to promote <a href="https://twitter.com/hashtag/reproducibility?src=hash&amp;ref_src=twsrc%5Etfw">#reproducibility</a> in <a href="https://twitter.com/hashtag/bioinformatics?src=hash&amp;ref_src=twsrc%5Etfw">#bioinformatics</a> analyses. <a href="https://twitter.com/alex_peltzer?ref_src=twsrc%5Etfw">@alex_peltzer</a> <a href="https://twitter.com/aka_hpatel?ref_src=twsrc%5Etfw">@aka_hpatel</a></p>&mdash; Rodrigo O. Polo (@ropolo) <a href="https://twitter.com/ropolo/status/1158860737516675073?ref_src=twsrc%5Etfw">August 6, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">Literally everyone should know about that! <a href="https://twitter.com/nextflowio?ref_src=twsrc%5Etfw">@nextflowio</a> Harshil Patel from Francis Crick Institute pitches good practice, hydrating and having these cool stickers on his laptop <a href="https://twitter.com/hashtag/bioinfocore?src=hash&amp;ref_src=twsrc%5Etfw">#bioinfocore</a> <a href="https://twitter.com/hashtag/ISMBECCB?src=hash&amp;ref_src=twsrc%5Etfw">#ISMBECCB</a> <a href="https://t.co/tkIgAWfeYZ">https://t.co/tkIgAWfeYZ</a> <a href="https://t.co/pOHNcFh3Yp">pic.twitter.com/pOHNcFh3Yp</a></p>&mdash; Daniel Stekhoven (@DanStekhoven) <a href="https://twitter.com/DanStekhoven/status/1153250437904130048?ref_src=twsrc%5Etfw">July 22, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">High throughput sequencing data processing keeps getting easier than ever thanks to <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a>! <a href="https://t.co/BsaYCOSq3I">https://t.co/BsaYCOSq3I</a></p>&mdash; Ignacio Tripodi (@ignaciot) <a href="https://twitter.com/ignaciot/status/1136688546251542528?ref_src=twsrc%5Etfw">June 6, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">Many useful and well-documented bioinformatics pipelines by the <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a> team, from RNA-seq, methylation, to this differential ATAC-seq. More in development, <a href="https://t.co/SvAJjIJqHD">https://t.co/SvAJjIJqHD</a> <a href="https://t.co/TCsBWlJyeC">https://t.co/TCsBWlJyeC</a></p>&mdash; Mikhail Dozmorov (@mikhaildozmorov) <a href="https://twitter.com/mikhaildozmorov/status/1116154350026489857?ref_src=twsrc%5Etfw">April 11, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">Alex from <a href="https://twitter.com/QBIC_tue?ref_src=twsrc%5Etfw">@QBIC_tue</a> talking about how <a href="https://twitter.com/nextflowio?ref_src=twsrc%5Etfw">@nextflowio</a> and <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a> enable his genomics pipelines, built on <a href="https://twitter.com/awscloud?ref_src=twsrc%5Etfw">@AWScloud</a> Batch. Really cool. <a href="https://t.co/FDuXdmem1n">pic.twitter.com/FDuXdmem1n</a></p>&mdash; Brendan Bouffler☁️🏳️‍🌈 (@boofla) <a href="https://twitter.com/boofla/status/1100764282617253888?ref_src=twsrc%5Etfw">February 27, 2019</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">My definition of joy: <br>1) Pick a random <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a> pipeline <br>2) Launch `nextflow run &lt;name&gt; -profile docker,test`<br>3) Feel proud of that.</p>&mdash; Paolo Di Tommaso (@PaoloDiTommaso) <a href="https://twitter.com/PaoloDiTommaso/status/1073568446238068736?ref_src=twsrc%5Etfw">December 14, 2018</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr">Just managed to get <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a> on <a href="https://twitter.com/nextflowio?ref_src=twsrc%5Etfw">@nextflowio</a> working on my machines. Really cool. Thanks for all your hard work</p>&mdash; Alastair Kerr (@alastair_kerr) <a href="https://twitter.com/alastair_kerr/status/1068491249315926016?ref_src=twsrc%5Etfw">November 30, 2018</a></blockquote></div>

<div class="card border-0"><blockquote class="twitter-tweet"><p lang="en" dir="ltr"><a href="https://twitter.com/tallphil?ref_src=twsrc%5Etfw">@tallphil</a> talking at the <a href="https://twitter.com/nextflow?ref_src=twsrc%5Etfw">@nextflow</a> delelopers’ hackathon about the incredible work he and his team have done to build <a href="https://twitter.com/nf_core?ref_src=twsrc%5Etfw">@nf_core</a>. <a href="https://t.co/SbWQIh2Q5r">pic.twitter.com/SbWQIh2Q5r</a></p>&mdash; Brendan Bouffler☁️🏳️‍🌈 (@boofla) <a href="https://twitter.com/boofla/status/1065614143061991424?ref_src=twsrc%5Etfw">November 22, 2018</a></blockquote></div>

</div>


<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


<?php include('../includes/footer.php');
