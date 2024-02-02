---
title: 'Bytesize 14: Graphic design / pipeline diagrams'
subtitle: Zandra Fagernäs - MPI-SHH
type: talk
startDate: '2021-05-25'
startTime: '13:00+02:00'
endDate: '2021-05-25'
endTime: '13:30+02:00'
youtubeEmbed: https://youtu.be/5jZPucWXnno
locationURL:
  - https://www.bilibili.com/video/BV1Z54y1V78h
  - https://youtu.be/5jZPucWXnno
  - https://doi.org/10.6084/m9.figshare.14673117.v1
---

This week, Zandra Fagernäs ([@ZandraFagernas](http://github.com/ZandraFagernas/)) will present: _**Graphic design / pipeline diagrams.**_

This will cover:

- How the graphics of nf-core/eager were created.
- Tips and tricks for how someone with no experience in illustration can make graphics.

The talk will be presented on Zoom and live-streamed on YouTube.

- YouTube: <https://youtu.be/5jZPucWXnno>
- Bilibili: <https://www.bilibili.com/video/BV1Z54y1V78h>

<details markdown="1"><summary>Video transcription</summary>
:::note
The content has been edited to make it reader-friendly
:::

[0:40](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=40) I would like to begin by thanking the organisers for inviting me to give this talk. I will be going through what we did for the nf-core/eager pipeline and why we made the specific choices that we did.

[1:01](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=64) I’m a PhD student at the Max Planck Institue for the Science of Human History (MPI-SHH), and I work with biomolecules in archaeological dental calculus. So I’m basically an end-user of the pipelines on nf-core. The reason I’m here is because I made the graphics for the `nf-core/eager` pipeline.

[1:30](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=90) But if I’m not an illustrator, then how am I qualified to do this? Well, I took a two-hour class in illustration around four years ago, which resulted in a [colouring book](http://christinawarinner.com/outreach/children/adventures-in-archaeological-science/).

[1:47](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=107) My mom and my grandma are really good illustrators. I’m not an artist at all, but maybe art is in my blood.

[1:57](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=117) I also have a Twitter page about what it’s like to work in an ancient biomolecules lab, where I have these illustrations that around seven other people find funny.

[2:10](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=130) I guess the main reason I got involved in the illustrations for eager was that I shared an office with James Fellows Yates (@jfy133) at the time, and James is one of the main driving forces behind `nf-core/eager`.

[2:26](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=146) So why did we choose to do the different graphics for `nf-core/eager`? Well I can tell you as an end user who is more comfortable in the lab, I think software documentation can be a little boring to read, especially if they are just long `README` files. I would likely just stop reading them after the second paragraph. These sorts of text files can also be very intimidating for a first-time user, and large pipelines like eager can be really hard to comprehend since they can have many different components.

[3:11](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=191) This is why the developers of `nf-core/eager` got me involved in creating illustrations for the documentation.

[3:26](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=205) James and I largely worked together. He would send me sketches like this one here, covering what he would like to have in the figure. Then I would spend around half an hour trying to decipher James’ handwriting, and looking at the actual figures that are outputs of `nf-core/eager`. I would then create a cutesy figure. I think our strength here was that James, who was very involved in writing the pipeline knew all the information that was required in the figure, and I was more focussed on illustrating those clearly.

[4:39](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=279) So for nf-core/eager, we decided to make two types of figures. The ones I will go through first serve to explain results, and the second ones are different sorts of pipeline diagrams.

[4:56](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=296) So for the figures in the documentation for `nf-core/eager`, I’m showing you an example of a figure that we have in our manuscript where panel “A” is the actual output from FastQC and panel “B” is the adaptation that we have made. We wanted to make them simple and easy to follow, so we removed a lot of information that was unnecessary for a representative illustration. We also wanted to make them more visually appealing through the use of slightly brighter colours, font types, and soft contours. But we also wanted to make them informative yet accessible, and encourage new users to try using the pipeline on their own. All these figures are available for use under a CC-BY licence, meaning that people can adapt them and use them for their own purposes as well.

[6:12](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=372) Then the second part is different pipeline diagrams, and what I’ve shown here is the main overview of the `nf-core/eager` pipeline. We wanted to make it simple and informative, and to be both self-explanatory and complement the documentation. So we have the different inputs marked by these different icons in square boxes to differentiate them from the pipeline elements, and the same with the outputs. We group the different elements of the pipelines in blocks in a way that is easy to follow for individuals running ancient DNA analysis. You can clearly see what happens after you provide your fastq file as an input; the different steps and choices you have at each step etc. The colours here indicate new additions to steps in the pipeline.

[7:23](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=443) I’d like to add here that I did not make these from scratch. James mentioned that he liked the pipeline diagram from Sarek, so I took a look and thought he was right. So we picked out elements that we liked, and this is what I suggest that everyone who wants to start making a pipeline diagram does. Just search for some on the internet, and remember that you don’t need to reinvent the wheel. You can just pick out things that you like and things you don’t from existing diagrams, and adapt them to make your own.

[8:02](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=482) I use Inkscape to make my diagrams, mainly because it is free and because it is intuitive. It doesn’t take long to learn to use and it has a strong online community that makes it easy to troubleshoot or ask for help. You can just search for solutions to a problem on the internet, and there are lots of online tutorials and DIY blog posts, which are helpful.

[8:47](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=527) I also wanted to mention colour schemes. I understand that these are often quite hard to pick. My suggestion is to keep it as simple as possible. I personally think that two-colour schemes are sufficient. For nf-core illustrations, just pick the codes from the logo, as you will likely be using the logo in your figure anyway. So if you pick the same colour scheme, you won’t be increasing the number of colours that you use.

[9:28](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=568) To illustrate the effects of using too many colours, take a look at this image here. I’ve used the colours from the `nf-core/eager` logo here, but I’ve also chosen to use blue to represent parts of the old pipeline, yellow to represent the different inputs, and the arrows in red. I think that makes it harder to focus and runs the risk of making it inaccessible to some users of the pipeline.

[9:58](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=598) If you’re wondering how to figure out if colours work well together, try and search for different colour schemes on the internet - there are tons of them - and pick out ones that are also accessible for users with colour-blindness.

[10:31](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=631) Finally, I want to mention that you are not alone in this. We have definitely relied a lot on the community to get feedback and inspiration for our figures. So at the [link](https://nf-co.re/docs/contributing/design_guidelines), you can find the nf-core logo, so you don’t need to cut it out from elsewhere, you can find the font-type that is used, and even the colours that are generally used. You can also find examples of previous pipelines and previous figures, icons for the different components such as fastq and psv icons that I used earlier. You can also always ask questions in `#graphics` on Slack. It is full of people who are happy to provide feedback on figures.

[11:19](https://youtu.be/5jZPucWXnno?list=PL3xpfTVZLcNiSvvPWORbO32S1WDJqKp1e&t=679) Thanks everyone. I have hopefully showed you that anyone can get started with illustrating and that one needn’t be a professional designer or an illustrator to make these sorts of diagrams.

</details>
