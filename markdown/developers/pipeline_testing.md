---
title: Pipeline testing
subtitle: Guidelines for testing.
---

# Testing with third-party tools

By default, every single pipeline has a `.travis.yml` file which defines testing
rules. This allows for consistency and reproducibility. Below, we will describe
platforms we use: Travis CI (cloud) and Jenkins (in-house only available for
SciLifeLab personnel).

# Travis CI

To be written by someone.

# Jenkins

We have build a custom **in-house cluster** consisting of one master and two
slaves. The master is responsible for running Jenkins dashboard while the slaves
are used for running the pipelines.

> NB: This cluster is only available to SciLifeLab users as it is behind
> firewall and reverse proxy.

## Installation

```bash
su
sh install.sh
docker-compose up -d
```

All installation scripts are included in [jenkins folder].

### Rebuilding stack

```bash
docker-compose down
docker-compose up --force-recreate --build -d
```

## Schema

- Master
  - Host: kraken
  - URL: [http://kraken.dyn.scilifelab.se](http://kraken.dyn.scilifelab.se)
  - Function: Jenkins master + reverse-proxy

- Nodes
  - Host: ship-1, ship-2
  - Function: Jenkins slaves

## Running services

### Traefik

Serves as load balancer/reverse proxy. In order to register new service, include
these labels in `docker-compose.yml`.

```yaml
labels:
  - traefik.port=8080
  - traefik.enable=true
  - traefik.backend=<SERVICE_NAME>
  - traefik.frontend.rule=Host:${http://kraken.dyn.scilifelab.se};PathPrefix:/<PREFIX>
  - traefik.frontend.passHostHeader=true
```

### Firewall (uwf)

```bash
ufw default deny incoming
ufw default allow outgoing
ufw allow ssh
# allow reverse-proxy
ufw allow 8080
ufw allow 80
# allow jenkins
ufw allow 50000
# Make sure to have an open SSH connection in case you have changed some rules
# otherwise you will lock yourself.
ufw enable
```

### Tips & tricks

1. Get admin password for Jenkins installation

`docker exec jenkins-master cat /var/jenkins_home/secrets/initialAdminPassword`

2. List all installed packages for automatic installation `jenkins/script`:

```groovy
Jenkins.instance.pluginManager.plugins.each{
  plugin ->
    println("${plugin.getShortName()}:${plugin.getVersion()}")
}
```

### Known issues

1. Jenkins blank page
Invoke restart by opening [restart link].

## Things to consider

- [Icinga2](https://icinga.com/docs/icinga2/latest/) (for monitoring)

[jenkins folder]: https://github.com/nf-core/nf-co.re/tree/master/includes/jenkins
[restart link]: http://kraken.dyn.scilifelab.se/jenkins/restart
