## Governance

nf-core is a community effort to collect a curated set of analysis pipelines built using Nextflow.

As a whole, nf-core is a community project. Anyone with an interest in the project can join the community and contribute.
To help manage the project, nf-core has teams that are made up of community members that oversee core activities.

<div class="mermaid">
flowchart BT
    core --> steering
    infrastructure --> core
    outreach --> core
    ambassadors --> outreach
    maintainers --> core
    safety --> steering
</div>

The purpose of this page is to formalize the governance for the nf-core project and community.

This page describes the governance structure of the nf-core community, how governance teams are elected, their responsibilities, and how decisions are made.

This page outlines the roles of the project teams, along with the responsibilities and privileges that come with them. It also outlines how team members are elected/chosen.

#### Steering committee {#steering}

Changes that impact the community require decisions informed by extensive experience with both the nf-core project and the larger ecosystem.
The steering committee includes representatives from the core team and outside advisors who oversee the running of the nf-core project. 

The steering committee is not a fixed size with its members being elected by the current steering committee based on need.
The steering committee will meet regularly to discuss the overall direction of the project, funds, and personnel.

**Responsibilities**

- Making decisions about project funds and personnel
- Guiding large initiatives

**Members**

<a class="btn btn-light rounded-pill"><i aria-hidden="true"></i>Ellen Sherwood</a>
<a class="btn btn-light rounded-pill" href="https://github.com/evanfloden"><i aria-hidden="true"></i>Evan Floden</a>
<a class="btn btn-light rounded-pill" href="https://github.com/ggabernet"><i aria-hidden="true"></i>Gisela Gabernet</a>
<a class="btn btn-light rounded-pill" href="https://github.com/drpatelh"><i aria-hidden="true"></i>Harshil Patel</a>
<a class="btn btn-light rounded-pill" href="https://github.com/ewels"><i aria-hidden="true"></i>Phil Ewels</a>

#### Core team {#core}

The core team ensures the day-to-day running of the nf-core project and oversees the activities of subgroups.

The core team is community members that have demonstrated a continued commitment to the nf-core community.
New members will be invited to be a part of the core team based on contributions, experience, and engagement with the community.
New members will be elected by current core team members.
The core team will aim to have representation from different genders, geography, and employment (e.g., academia, clinical and industry).

Significant community decisions will be made by vote with any decision without a clear majority being passed to the steering committee.

Core team members will appear as organization members on the GitHub organization and have administrator access to repositories.

**Responsibilities**

- Day-to-day community decisions
- Attendance at the core team annual meeting
- Regular attendance at monthly core team meetings
- Sub-roles within the nf-core governance structure
- A strong community presence

**Members**

<a class="btn btn-light rounded-pill" href="https://github.com/apeltzer"><i aria-hidden="true"></i>Alexander Peltzer</a>
<a class="btn btn-light rounded-pill" href="https://github.com/christopher-hakkaart"><i aria-hidden="true"></i>Chris Hakkaart</a>
<a class="btn btn-light rounded-pill" href="https://github.com/FriederikeHanssen"><i aria-hidden="true"></i>Friederike Hanssen</a>
<a class="btn btn-light rounded-pill" href="https://github.com/ggabernet"><i aria-hidden="true"></i>Gisela Gabernet</a>
<a class="btn btn-light rounded-pill" href="https://github.com/drpatelh"><i aria-hidden="true"></i>Harshil Patel</a>
<a class="btn btn-light rounded-pill" href="https://github.com/jfy133"><i aria-hidden="true"></i>James A. Fellows Yates</a>
<a class="btn btn-light rounded-pill" href="https://github.com/JoseEspinosa"><i aria-hidden="true"></i>Jose Espinosa-Carrasco</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mirpedrol"><i aria-hidden="true"></i>Júlia Mir Pedrol</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mribeirodantas"><i aria-hidden="true"></i>Marcel Ribeiro-Dantas</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mashehu"><i aria-hidden="true"></i>Matthias Hörtenhuber</a>
<a class="btn btn-light rounded-pill" href="https://github.com/MaxUlysse"><i aria-hidden="true"></i>Maxime Garcia</a>
<a class="btn btn-light rounded-pill" href="https://github.com/ewels"><i aria-hidden="true"></i>Phil Ewels</a>

**Core team alumni**

<a class="btn btn-light rounded-pill" href="https://github.com/alneberg"><i aria-hidden="true"></i>Johannes Alneberg </a>
<a class="btn btn-light rounded-pill" href="https://github.com/KevinMenden"><i aria-hidden="true"></i>Kevin Menden</a>
<a class="btn btn-light rounded-pill" href="https://github.com/olgabot"><i aria-hidden="true"></i>Olga Botvinnik</a>
<a class="btn btn-light rounded-pill" href="https://github.com/renbot-bio"><i aria-hidden="true"></i>Renuka Kudva</a>
<a class="btn btn-light rounded-pill" href="https://github.com/sven1103"><i aria-hidden="true"></i>Sven F.</a>

#### Safety {#safety}

The nf-core community should feel comfortable contributing to the project without the risk of harassment or abuse.
The safety team is responsible for ensuring the community is a safe place and responding to instances of misconduct.

The safety team is community members who have displayed integrity, strong communication, and a genuine concern for community welfare and are selected by the core team.
The safety team is not a fixed size and will scale as the community grows.
The safety team is not a part of the core team and can report directly to the steering committee.

**Members**

<a class="btn btn-light rounded-pill" href="https://github.com/ctuni"><i aria-hidden="true"></i>Cris Tuñí</a>
<a class="btn btn-light rounded-pill" href="https://github.com/heuermh"><i aria-hidden="true"></i>Michael Heuer </a>
<a class="btn btn-light rounded-pill" href="https://github.com/snafees"><i aria-hidden="true"></i>Saba Nafees</a>

**Responsibilities**

- Write and maintain the nf-core code of conduct
- Be available for nf-core events (online or in person)
- Monitor Slack channels for instances of misconduct
- Promptly respond to reports of misconduct and escalate to the core team or steering committee as necessary

#### Infrastructure {#infrastructure}

Tooling is a fundamental part of the nf-core community.
The infrastructure team is responsible for the development and implementation of the nf-core tooling framework.
The infrastructure team will have one or more leads who are responsible for overseeing infrastructure efforts.
Infrastructure team members are elected by the core team.
The infrastructure team is not a fixed size and will scale as the community grows.

**Responsibilities**

- Development and maintenance of nf-core tools
- Development and maintenance of the nf-core website
- Development and maintenance of nf-core megatests
- Regular attendance at maintenance team meetings

**Members**

<a class="btn btn-light rounded-pill" href="https://github.com/ewels"><i aria-hidden="true"></i>Phil Ewels</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mirpedrol"><i aria-hidden="true"></i>Júlia Mir Pedrol</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mashehu"><i aria-hidden="true"></i>Matthias Hörtenhuber</a>

#### Outreach {#outreach}

Outreach is an important part of any community project.  
The outreach team is responsible for overseeing the organization and running community outreach efforts, including, but not limited to, the nf-core ambassador program, hackathons, and the `#bytesize` seminar series.

The outreach team will have one or more leads who are responsible for overseeing outreach efforts.
New members will be invited to be a part of the outreach team based on experience and outreach activity (e.g., involvement in the ambassador program). 
The outreach team is not a fixed size and will scale as the community grows.

The outreach leads will have access to community social media and YouTube accounts (e.g., Twitter and YouTube).

**Responsibilities**

- Organizing and running the `#bytesize` seminar series
- Leading the organizing hackathons and training
- Organizing the ambassador program
- Creating and sharing community content
- Regular attendance at outreach team meetings

**Members**

<a class="btn btn-light rounded-pill" href="https://github.com/christopher-hakkaart"><i aria-hidden="true"></i>Chris Hakkaart</a>
<a class="btn btn-light rounded-pill" href="https://github.com/FranBonath"><i aria-hidden="true"></i>Franziska Bonath</a>
<a class="btn btn-light rounded-pill" href="https://github.com/mribeirodantas"><i aria-hidden="true"></i>Marcel Ribeiro-Dantas</a>
<a class="btn btn-light rounded-pill" href="https://github.com/abhi18av"><i aria-hidden="true"></i>Abhinav Sharma</a>
<a class="btn btn-light rounded-pill" href="https://github.com/Miller"><i aria-hidden="true"></i>Edmund Miller</a>
<a class="btn btn-light rounded-pill" href="https://github.com/MaxUlysse"><i aria-hidden="true"></i>Maxime Garcia</a>
<a class="btn btn-light rounded-pill" href="https://github.com/pcantalupo"><i aria-hidden="true"></i>Paul Cantalupo</a>
<a class="btn btn-light rounded-pill" href="https://github.com/yuukiiwa"><i aria-hidden="true"></i>Yuk Kei</a>

#### Maintainers {#maintainers}

nf-core test data, modules, and pipeline repositories require regular upkeep and maintenance.
The maintainer's team is responsible for the management of these repositories in collaboration with the wider nf-core community.
The maintainer's team will have one or more leads who are responsible for overseeing maintenance efforts.
New members are invited to be a part maintainers team by current maintainers based on experience and community activity. The maintainer's team is not a fixed size and will scale as the community grows.

nf-core maintainers will have write access to repositories.

**Responsibilities**

- Respond to `#github-invitations`
- Review module and subworkflow pull requests
- Review pipeline releases (including first releases)
- Manage repository access for community developers
- Manage test data
- Enable and promote Flux community values

**Members**

<a class="btn btn-light rounded-pill" href="https://github.com/apeltzer"><i aria-hidden="true"></i>Alexander Peltzer</a>
<a class="btn btn-light rounded-pill" href="https://github.com/FriederikeHanssen"><i aria-hidden="true"></i>Friederike Hanssen</a>
<a class="btn btn-light rounded-pill" href="https://github.com/ggabernet"><i aria-hidden="true"></i>Gisela Gabernet</a>
<a class="btn btn-light rounded-pill" href="https://github.com/drpatelh"><i aria-hidden="true"></i>Harshil Patel</a>
<a class="btn btn-light rounded-pill" href="https://github.com/jfy133"><i aria-hidden="true"></i>James A. Fellows Yates</a>
<a class="btn btn-light rounded-pill" href="https://github.com/JoseEspinosa"><i aria-hidden="true"></i>Jose Espinosa-Carrasco</a>
<a class="btn btn-light rounded-pill" href="https://github.com/MaxUlysse"><i aria-hidden="true"></i>Maxime Garcia</a>

<a class="btn btn-light rounded-pill" href="https://github.com/Emiller88"><i aria-hidden="true"></i>Edmund Miller</a>

#### Ambassadors {#ambassadors}
