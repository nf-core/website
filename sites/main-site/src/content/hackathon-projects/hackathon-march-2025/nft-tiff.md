---
title: Continuation of nft-tiff
category: components
intro_video: ""
slack: https://nfcore.slack.com/archives/C07PTBHLK7D
image: "/assets/images/events/2025/hackathon-march/nft-tiff_meme.jpg"
image_alt: Remote sensing pipelines without nft-tiff are a pain to test!
leaders:
  Felix-Kummer:
    name: Felix Kummer
    slack: "https://nfcore.slack.com/team/U04SDFLTBA5"
---

This project aims to continue the development of the nf-test plugin [nft-tiff](https://github.com/nf-core/nft-tiff) for `.tiff`-files.

## Goals

1. Implement exact matching for `.tiff` files (content and metadata).
2. Implement partial matching for `.tiff` files (e.g. x % agreement).

## Motivation

The need for [nft-tiff](https://github.com/nf-core/nft-tiff) arose during the development of [nf-core/rangeland](https://nf-co.re/rangeland), a pipeline that analyses satellite imagery.
As the first earth observation pipeline in nf-core with a stable release, no off-the-shelf nf-test plugin for checking common earth observation file formats was available.
The development of nft-tiff started at a previous hackathon.
For this hackathon, we would like to bring nft-tiff closer to a full release.
The development of this plugin will benefit future earth observation pipelines in nf-core and Nextflow as it supports a common file format in the field.

## Tagged Image File Format (TIFF) and GeoTIFF

TIFF is a file format for storing multi-channel raster images.
An important extension to the TIFF format is GeoTIFF.
GeoTIFF extends the format with a variety of metadata fields relevant for geoscientific analysis (such as [nf-core/rangeland](https://nf-co.re/rangeland)).
