<?php
$title = 'About nf-core';
$subtitle = 'Details about the nf-core project - who is involved and how it was started.';
include('../includes/header.php');
?>

<h1 id="about-nf-core">About nf-core</h1>
<p>The nf-core project came about at the start of 2018 as the pet-project of <a href="http://phil.ewels.co.uk/">Phil Ewels</a>. He’s the head of the development facility at <a href="https://ngisweden.scilifelab.se/">NGI Stockholm</a> (National Genomics Infrastructure), part of <a href="https://www.scilifelab.se/">SciLifeLab</a> in Sweden.</p>

<p>The NGI has been developing analysis pipelines for use with it’s genomics data for several years and started using a set of standards for each pipeline created. This helped other people run the pipelines on their own systems - typically Swedish research groups at first, but later on other core genomics facilities too.</p>

<p>As the number of people interested in these common pipelines grew, it seemed silly to keep all the pipelines under the umbrella of just SciLifeLab, complete with NGI prefixes and logos. To try to open up the effort into a truly collaborative project, <a href="https://github.com/nf-core">nf-core</a> was created and all relevant pipelines moved to this new GitHub Organisation.</p>

<h2 id="participating-groups">Participating Groups</h2>
<p>A number of groups are involved in running, maintaining and developing pipelines that are part of nf-core. They are:</p>

<ul>
  <li>The <a href="https://ngisweden.scilifelab.se/">National Genomics Infrastructure</a> (NGI), SciLifeLab Sweden</li>
  <li>The <a href="qbic.life">Quantitative Biology Center</a> (QBiC), Universität Tübingen, Germany</li>
  <li>The <a href="http://www.crg.eu/">Centre for Genomic Regulation</a> (CRG), Barcelona, Spain</li>
  <li>The <a href="https://www.a-star.edu.sg/gis">Genomics Institute of Singapore</a> (GIS), A*STAR, Singapore</li>
  <li>The Competence Centre for Genome Analysis Kiel, Kiel University, Germany</li>
  <li>The <a href="https://www.iarc.fr">International Agency for Research on Cancer</a> (IARC), World Health Organisation, Lyon, France</li>
</ul>

<p>Is your group missing? Please <a href="https://github.com/nf-core/nf-core.github.io/blob/master/about.md">submit a pull request</a> to add yourself!</p>

<h1 id="contributors">Contributors</h1>
<p>nf-core is by design a collaborative effort, and would not exist if it were not for the efforts of many dedicated contributors. They include:</p>

<ul>
  <li>Phil Ewels (<a href="https://github.com/ewels/">@ewels</a>) - Project lead</li>
  <li>CRG Barcelona
    <ul>
      <li>Paolo Di Tommaso (<a href="https://github.com/pditommaso">@pditommaso</a>) - Developer of Nextflow and source of never-ending support</li>
      <li>Evan Floden (<a href="https://github.com/skptic">@skptic</a>)</li>
    </ul>
  </li>
  <li>NGI Sweden
    <ul>
      <li>Rickard Hammarén (<a href="https://github.com/hammarn">@Hammarn</a>)</li>
      <li>Chuan Wang (<a href="https://github.com/chuan-wang">@chuan-wang</a>)</li>
      <li>Denis Moreno (<a href="https://github.com/Galithil">@Galithil</a>)</li>
      <li>Remi-Andre Olsen (<a href="https://github.com/remiolsen">@remiolsen</a>)</li>
      <li>Maxime Garcia (<a href="https://github.com/MaxUlysse">@MaxUlysse</a>)</li>
      <li>Senthil (<a href="https://github.com/senthil10/">@senthil10</a>)</li>
    </ul>
  </li>
  <li>QBiC Tübingen
    <ul>
      <li>Alexander Peltzer (<a href="https://github.com/apeltzer">@apeltzer</a>)</li>
      <li>Sven Fillinger (<a href="https://github.com/sven1103">@sven1103</a>)</li>
    </ul>
  </li>
  <li>GIS Singapore
    <ul>
      <li>Andreas Wilm (<a href="https://github.com/andreas-wilm">@andreas-wilm</a>)</li>
      <li>Chih Chuan (<a href="https://github.com/shihcc">@shihcc</a>)</li>
    </ul>
  </li>
  <li>Competence Centre for Genome Analysis Kiel
    <ul>
      <li>Mark Hoeppner (<a href="https://github.com/marchoeppner">@marchoeppner</a>)</li>
    </ul>
  </li>
  <li>International Agency for Research on Cancer
    <ul>
      <li>Matthieu Foll (<a href="https://github.com/mfoll">@mfoll</a>)</li>
      <li>Tiffany Delhomme (<a href="https://github.com/tdelhomme">@tdelhomme</a>)</li>
    </ul>
  </li>
</ul>

<p>Want to be added to this list? Please <a href="https://github.com/nf-core/nf-core.github.io/blob/master/about.md">submit a pull request</a> to add yourself!</p>



<?php include('../includes/footer.php'); ?>
