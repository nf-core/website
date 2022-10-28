---
title: 'Bytesize: nf-core/inkscape and tube map diagrams'
subtitle: James Fellows Yates - Max Planck Institute for Evolutionary Anthropology, Leipzig, Germany
type: talk
start_date: '2022-07-12'
start_time: '13:00 CEST'
end_date: '2022-07-12'
end_time: '13:45 CEST'
youtube_embed: https://www.youtube.com/watch?v=0vKhfedYKGo
location_url:
  - https://www.youtube.com/watch?v=0vKhfedYKGo
  - https://doi.org/10.6084/m9.figshare.20296869.v1
---

# nf-core/bytesize

Join us for our **weekly series** of short talks: **“nf-core/bytesize”**.

Just **15 minutes** + questions, we will be focussing on topics about using and developing nf-core pipelines.
These will be recorded and made available at <https://nf-co.re>
It is our hope that these talks / videos will build an archive of training material that can complement our documentation. Got an idea for a talk? Let us know on the [`#bytesize`](https://nfcore.slack.com/channels/bytesize) Slack channel!

## Bytesize: nf-core/inkscape and tube map diagrams

This week, James Fellows Yates ([@jfy133](https://github.com/jfy133)) will give a basic introduction to the vector-based graphic design tool [Inkscape](https://inkscape.org/), and then how to use the tool to create the popular 'tube map'-style pipeline diagrams.

![nf-core/eager pipeline diagram in the style of a tube-map](https://raw.githubusercontent.com/nf-core/eager/master/docs/images/usage/eager2_metromap_complex.png)

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
(host) Hello everyone, I'm here today with James Fellows Yates and he's going to talk about Inkscape and pipeline diagrams.

Thank you very much. This is a repeat of a talk I gave back at the last nf-core hackathon and what I'll be talking about is how you can make pretty diagrams aimed at improving your documentation. This can be making it more eye-catchy, more accessible and making people understand better how to run your pipeline and how to interpret it. What I'll do is give a very brief introduction to what things you can generate with a tool called Inscape. Then I'll give a live demo of how you use Inkscape and the main components and functionality that I use to make and the different pipeline map constructions, which become semi-famous, at least on Twitter. Then I'll show you how I will would construct a pipeline diagram such as this.

[1:08](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
The tool that we're talking about today is Inkscape, which is an open source, pretty close to professional quality software, which can run on Linux OSx and Windows, which is really nice. It's equivalent to to Adobe Illustrator, I think, but it is but open source, so it's free and accessible to everybody. It's a really nice and quite polished tool and I really enjoy it a lot and I've been using it for quite a few years now. The reason why we want to use Inkscape is primarily, because of the difference between raster and vector images. When you take a photo on your phone or with a digital camera, the images that you generate are raster images, so if you zoom in really closely, you start getting these very pixely squares of different colors. 

[1:58](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
This is very nice in many ways, but when you're dealing with the web, where a lot of bioinformatic documentation is displayed on (and we people like to zoom in and out and stretch and modify things) and to make sure like your images look very nice, clean and clear... When you're dealing with particular web documentation, we prefer a format called vector images. In particular the one we'll be using today is SVG. The really nice thing is that you can zoom, stretch, manipulate and you always have these smooth, very high definition lines. The way these SVGs are constructed, which are based on coordinate system, is that they are very portable. You can take different components, break it up of a particular image, and import it into other images. Examples of vector images... Things you can do to improve your documentation in terms of graphic design, when it comes to bioinformatic pipelines like in nf-core, is for example to design a logo. Here's a couple of examples of nf-core: there is Sarek and Eager is another one, that I made for a different project. Logos are a very good way of making your pipelines look a bit more professional and more eye-catching. People are more likely going to trust the quality of the pipeline when it has a brand identity, which is a nice little productive procrastination if you ever need that.

[3:34](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
You can also make workflow diagrams. This is one from Eager. This is a more a general, broad overview, which helps people to understand what are the major components of the pipeline and what it does. It doesn't go into too much detail and these are the sections and the tools it uses. If you want something more detailed, you can make something like this, which shows more the direct connections between the different components, the different tools of the particular pipeline. These give you two different purposes: This is more for like selling the pipeline, if you want to pique someone's interest, by for example putting by it in a publication. Whereas this is more for the user who's actually going to be interacting with the pipeline, trying to run it. There are different stages and levels of such pipeline diagrams. You can use them to help make your pipeline more accessible from the concept of what it does.
This format is the one that became quite popular recently, but i've been playing around with other variants. This is a new one, that I'm coming up with for nf-core text profiler. This shows that you don't have to stick with the particular one that I showed in a previous slide, you can play around with this as well.

[4:51](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
All of these diagrams, all these images are done in Inkscape. In addition, it doesn't have to just be diagrams to show what you can do in your documentation. One thing that we did with Eager, which was also quite popular and I was quite proud of, is that my co-author Sandra Wegener, she drew these very nice schematics cartoony diagrams for the output documentation. This helps people to understand better what they should be looking for when they're doing quality control in their pipeline runs. For example, she took what you get from MultiQC, so for example FastQC in this case, the first QC section. She drew little cartoony versions of what you should be seeing, separated by the three different boxes, but with little notes saying exactly what you are seeing in this and how you should interpret the results here. This acts as a very nice, quick and fast reference for people, to understand what your pipeline has done at the other end. Also it adds a little bit of fun, a little bit of color, to otherwise often very dry documentation, that can be often be very dull. Particularly for very big pipelines like Eager and Sarek, and when you have a lot of output documentation, such things can break things up a little bit and make it more accessible and more fun for a user to use.

[6:20](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
There are other things you can do with Inkscape. For example these are a couple of other things I've done. This was a a schematic overview of some characteristics of ancient DNA, done all in Inkscape. Actually all of these little Emojis, they have been imported from another project with OpenMoji. Everything is in SVG and it's very easy to drag and drop it into my image, recolor things, change the size and so on, which is really nice. That's all benefits from this vector imaging format! Also you can do more realistic drawings by tracing and, doing this in Inkscape. It's not so cartoony as all the other previous objects. This is a tooth, cut in half, and that's showing how to sample for ancient DNA. You can load the raster image into Inkscape then trace over the top and then fill in the colors. As such you can get a particularly good representation of the object, which is then very manipulable after.

[7:23](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
That's the things you can do with Inkscape. Now, what i'd like to do, is to go into Inkscape itself and show you some of the basic functionalities that you would need, in order to make one of the pipeline diagrams, such as this one and this one.

[7:44](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
When you open Inkscape it should look something like this.
The UIs (user interface) should be relatively familiar: you've got the toolbar right at the top, you've got other toolbars on the sides as well. It may look a little bit different, depending on what operating system you are you're working on, but generally this is the outline you will have, at least in the latest version. What you should hopefully see on my screen is a a key logger, you should be able to see what I'm typing. In the case I use like a shortcut, you can also see that, so keep an eye on that as well. 

[8:23](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
The most basic thing you can do in Inkscape, or with a vector image program, is making objects and shapes. I'm going to make here an object, which is a circle. I can make two circles, I can make squares, you can make triangles and convert these to stars and things like hexagons, or whatever. You can also resize these. As you can see, when I click on an object, it often will have either arrows like this, which you can resize like so.Sometimes if you click on this path editor (which I'll explain in a second) you can also edit with these squares and triangles here.

[9:05](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
You don't have to make fixed shapes like that. You can do lines with a the Bezier tool here. If I click here and click here, I'll make a line there; and with this Bezier tool, you're generating things called paths. These objects are fancy paths, but you can also create your own, random shapes you would like to make. For example like this. You're also not fixed to having sharp corners or sharp edges. You can bend edges and you can modify the corners with the options up here, to make them curved as well. Often this path system is how you do a tracing map, as I explained in the previous slide. You can also modify these as well. You can break them apart at the nodes for example... 
... and then join them together again. There's many many different ways you can manipulate these. 
I should clear up slightly... 

[10:21](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
Another thing you can do is to group objects together. You have these two objects and, let's say, you're happy with the way they are positioned now, and you want to keep this relationship as you manipulate the image. We can click both: holding shift and clicking the second one. You can press CTRL + g to move it together. You can see now, I'm moving both at the same time, but they're not losing their position. Furthermore, all objects have an order, they're overlapping each other. This order you can modify. 
That's probably a bad color pairing, sorry, one second, I'll just try this. 
What you can see here is, that the two rings overlap each other. By pressing page up and page down, I'm moving which is displaying on top. These buttons are here as well, which does the same thing. You can also move to the very top of the stack so you can have as many objects as you want. For example, I can send this pink one right to the bottom as well, like this.

[11:30](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
Another thing you can do, because everything's based on a coordinate system, is semi-programmatically reorder everything: aligning and distributing them based on some rule. Let's say I select all of these objects and I want to put them all in a single line, if I go to `Object` and then `Align and Distribute`, which should open a panel on the right hand side here. Your have various options on how to arrange and distribute your different objects. Let's say I want to have them all in a horizontal line, I can press this button here `Center on Horizontal Axis` and they'll all go into the line here. The same thing also goes for a vertical distribution.

[12:08](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
Maybe you want to have everything somewhat equally spaced out. You can see these are closer than these ones, and so on. You can do this with these buttons here on `Distribute`. Then I can blend it, like that, and you can see now, that everything is pretty much ordered in exactly the same way.

[12:29](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
These objects you can also scale and transform. Scaling. Like I said before, you can drag with the arrows here, but you can also go to the `Transform` under the `Object` menu here.  You can scale it up programmatically again. I want this one to be twice the size, and press apply, you see twice the size now. If I get the square... You can also rotate things as well, so let's say 45 degree angle, to get a sort of diamond

[13:00](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and you can also do text and all the text generally is also vector based so you can re reorder it detach them separate them out and so on so again clicking this text tool here I can start writing something I will say that the latest version of Inkscape because I think i've actually accidentally saw the beta version is a bit slow with the with the text for some reason I'm not entirely sure why so it might be a bit slower here but you get the idea you can make it bigger like this and and you can also color in the same way as the objects as well

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
so another nice thing an important thing for the pipeline diagrams is you often want to be able to as well as align and distribute things as you were doing here is bind them together in a consistent object together in a consistent manner so if you press the hash key on your keyboard you can get a grid and also in the settings under docking properties you can change the size of the grid and the style and with this it makes it easier for you to position things correctly and also when you have a grid on you can turn on a snapping tool as you can see up up here that you can also turn on with the percent key and with this you can have different types of edge snapping so what I mean by that is let's say I have turned on the midpoint snapping you can see now this little square coming popping up and is basically moving your object or the middle of the object in this case to the closest point on the grid here the same goes for the in the top corner here you have to get a feeling for exactly where to hold the object but basically you'll snap it so this way I can make sure for example manually I can see on the grid I want to have everything along this line here and obviously use the snapping to basically alter that accordingly

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and with coloring so you can see here in these objects currently I have only cut the outline of the object this is called a stroke you can also fill in as well so for example if I click oh no sorry

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
if I go to file then fill in stroke uses another menu here and you can set the colors or not why is it not working because there we go they're transparent so all objects can be transparent as well

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
you can see that basically there's two colors the outer line is called stroke and inside the object you can have the fill and these can be independently set which is very nice and this applies to any object want you can change the colors based on hex and colors if you want rgb wheels you can also set different stroke styles so at the moment let's change the color again something more obvious you can use a solid line you can also change the thickness of the of the outline and you can also change the style in terms of dashes let's move that example like so so these are all lots of flexibility in this in this manner as well so for example if you wanted to go with a web-based color you can like say you can use these hex im hex codes down the bottom so as an example if I were to go to the nfl graphics guidelines page just under the documentation and contributing you can see for example here under fonts and colors we have the official rgb values and also the hex value for this particular the end of call green so if I switch back here I can change the stroke of this object to the green as so again or rather the film in this case sorry oop and turning off the transparency and now this is the nf-4 green sometimes in some cases you may not want an svg as the final image you want to have in a different format that is also possible to to do so if you go to file and export you have the option to export for example as a png image you can just change it at different resolutions depending on whether you're going to be printing or put a web and things like that as well and this means that you also have that flexibility for different formats

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and how I'm doing for time okay now we can move on to the question which is how would you actually create the pipeline diagrams now you could create all the different components yourself separately but that speaking from experience the first time I did this that can take a lot of time so what we made for the mfcor website on the game the graphics design guidelines is actually a cheat sheet with all the different components you can use to make a pretty even even evenly spaced and distributed pipeline diagram we also have other components here you can download for example example of other pipelines if you want to modify these you can check the licenser over here there's different components for example for saric and maxine made very nice file icons which you can also download and using your own pipeline diagrams but for this example I'll just use the ones for the pipeline components which I made before so I can simply save this as a file

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and then drag this particular file into my Inkscape here and press okay and you can see here now i've got all of the objects there so like I emphasized before the nice thing about svg images because they're basically based on coordinate systems they can be rendered in many different ways and also downloading and importing this file into my document means I can actually now detach the different components of this and reuse them so for example I could take let's say my starting point and drag this over here if I follow the grid I will snap it with the midpoint snap here then I could take the straight line here and drag this over here and stick it on there I could put a another station and put that here but let's say now I want to do a split so I want to have basically two different lines because i've let's say optional pipelines or mutually exclusive optional steps to the pipeline I can copy and paste this over here and then by pressing h or v I can also flip and rotate the objects in that manner so I will again snap with the midpoint here and the same thing here and connect that there

oops

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
you get the point an important thing to to emphasize as you may have already noticed I'm just copying pasting and like i'd say 80 of my Inkscape work is basically copying and pasting components i've ever already made so that is why we made this cheat sheet so I highly recommend you basically be using this as well but of course you can modify them in whatever way you want but to get a basic outline i'd recommend trying something like this and so modifying let's say I don't like this green color I prefer to ungroup this I want a different color I can say change this to a red I'm gonna appear I would like this to have a purple on this one to be a blue blue and basically just keep working on this and construct things in this very lego-like fashion is my main recommendation in this case I can then also basically add a station name I'll call it let's say input put this up here put this down here say step one you want to change the fonts of your of your text so for example if you want to use the nf core font we have the things called maven pro exactly you can set this here as well

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and the final thing is maybe once you've completed your your image you want to not have this weird a4 size page you can also actually go to file looking properties and modify the layout in this way so for example setting is horizontal sorry landscape and also you actually resize the page to content if I press this button here

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
you can see that basically the page has been resized to cover all the different drawings that you may have in the particular document

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
and this is pretty much it and the one other thing I would recommend doing before actually starting such a diagram is already have a working doodle or diagram of your of your pipeline because sometimes this can get a bit fiddly when you're trying to work out exactly where everything should be spaced so for example for nf tax profiler we used google drawings and to basically sketch out already all of the different parts of the pipeline it's a bit easier to move things around here by snapping them together and then you can use this basically as a reference to how to make your your actual final work for documentation for example you could take the example at the top if I can go I said is the input I have the two splits here then i'd set up like this and basically follow along the guidelines so I do that here

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
so that is pretty much it that is the basics of what you need to to learn and so to recap you have all these objects you can make but generally when you come to the pipeline diagrams it's better to try and reuse components already made copying and pasting is your friend and you can color things by the menu and objects and fill in stroke and you have the fill which is the the outline of the object and the sorry the stroke of the outline the object and the fill is the inside of the object you can stretch manipulate all of the objects either manually by these handles or also programmatically with the transform and line and distribute objects so otherwise also have fun doing it again it's productive procrastination it does make a influence a lot of people call their eye for example the eager pipeline and because of such workflow diagrams we have a graphics design channel on slack which you also call these graphics which you can join as well if you have any questions but otherwise I think that is it so thank you very much are there any questions thank you very much so are there any questions in the audience

[0:01](https://www.youtube.com/watch?v=0vKhfedYKGo&t=1)
i have also allowed yourself to allow to unmute yourself now if anyone wants to say something otherwise thank you very much james and I want to also thank the john zuckerberg initiative for funding these talks and as usual you can continue the discussion on bite size slack channel bite size yes and thank you again