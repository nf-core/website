---
title: "Bytesize: nf-core/proteinfamilies"
subtitle: Evangelos Karatzas, EMBL Cambridge
type: talk
startDate: "2025-04-08"
startTime: "13:00+02:00"
endDate: "2025-04-08"
endTime: "13:30+02:00"
locations:
  - name: Online
    links:
      - https://kth-se.zoom.us/j/68390542812
---

In this weeks bytesize talk, Evanglos ([@vagkaratzas](https://github.com/vagkaratzas)) is going to introduce the nf-core pipeline nf-core/proteinfamilies.

The nf-core/proteinfamilies pipeline generates protein families from amino acid sequences and/or updates existing families with new sequences.
It takes a protein fasta file as input, clusters the sequences and then generates protein family Hiden Markov Models (HMMs) along with their multiple sequence alignments (MSAs).
Optionally, paths to existing family HMMs and MSAs can be given in order to update with new sequences in case of matching hits.
