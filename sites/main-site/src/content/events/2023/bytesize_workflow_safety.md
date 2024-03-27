---
title: 'Bytesize: Workflow safety and immutable objects'
subtitle: Rob Syme, Seqera Labs
type: talk
startDate: '2023-05-23'
startTime: '13:00+01:00'
endDate: '2023-05-23'
endTime: '13:30+01:00'
youtubeEmbed: https://www.youtube.com/watch?v=A357C-ux6Dw
locations:
  - name: Online
    links: https://www.youtube.com/watch?v=A357C-ux6Dw
---

This week, Rob Syme ([@robsyme](https://github.com/robsyme)) will talk about how to avoid introducing subtle concurrency bugs in your Nextflow workflows by safely modifying objects (or more specifically, by _not_ modifying objects) when writing Groovy closures.

<details markdown="1"><summary>Video transcription</summary>
**Note: The content has been edited for reader-friendliness**

[0:01](https://www.youtube.com/watch?v=A357C-ux6Dw&t=1)
(host) Hello everyone, Maxime here. First of all, I'd like to thank the Chan Zuckerberg Initiative to help us doing these bytesize talks. Today, I'm introducing Rob Syme, a Groovy wizard, working at Seqera Labs, and he will explain to us how to safely handle metadata. As usual, you'll be able to unmute yourself at the end of the talk to ask your question. Of course, you can also ask them on Slack or on the chat. Over to you, Rob.

[0:39](https://www.youtube.com/watch?v=A357C-ux6Dw&t=39)
Perfect. Thanks for the intro, Maxime. This talk is, as Maxime indicated, about safe metadata handling in Nextflow pipelines. But really, it's all about safely passing any objects through channels in Nextflow. But it has particular application for the nf-core community here because we do a lot of metadata handling in nf-core, particularly passing around meta-maps between processes and through channels. There are some complicated little bugs that can occur when mutating those objects in place. What I'm going to suggest to you is that you should never mutate those objects in place, but instead always return new objects. Of course, this is a bit hard to explain in the abstract. Let's do an example.

[1:23](https://www.youtube.com/watch?v=A357C-ux6Dw&t=83)
This example is based on a true story, a sad story. Let's take this workflow. We have a channel that produces an object or emits an object, and we're passing that channel to a BuyNewShirt process. This process does something very simple. It just echoes the t-shirt size and based on the weight property in that object being passed to the process. When we run this, we see the t-shirt size is small because the weight is 70 here, pass it to BuyNewShirt. That's what we'd expect. No surprises here. Nothing outrageously complicated. These weights are in kilos still. Let's make this a little bit more complicated. I'm assuming that everyone has a pretty solid grasp of this. Let's add two new processes. The one we're going to be specifically interested in here is GetNewJob. This GetNewJob process takes an object, emits that same object, and modifies that object in place. It adds five to the weight property. We take that channel again, the same channel that we created before, we pass it once to BuyNewShirt and once to GetNewJob. This runs as expected. The t-shirt size is small. Everything is as we expect. I think the way that a lot of people think about these pipelines is like this. I think conceptually a lot of us have this picture in our head when we're running these pipelines. We take an object, we pass it to BuyNewShirt, and we pass it to GetNewJob. We pass the object, those two channels into new processes. I think we have a tendency to think about those as independent events, but really because this channel is taking the same object, it's exactly the same object in memory getting passed to GetNewJob and BuyNewShirt. If GetNewJob modifies the properties of that, it will affect BuyNewShirt. In this case, we saw it didn't affect the output of the run because this is dependent on timing. It's an important understanding to remember that all of these processes happen asynchronously in Nextflow. We can offer no guarantees about the order or the timing of this.

[3:54](https://www.youtube.com/watch?v=A357C-ux6Dw&t=234)
To make it a little bit more clear, let's add in a little delay. I'm now going to add this Browse process before BuyNewShirt. This Browse process, all it does is sleep for a couple of seconds. That gives it time for the sort of... just delays before we BuyNewShirt. Now we've just added a process here. Ideally I think a lot of people conceptually would think that adding a process here shouldn't affect the outcome here because we've just passed the object through transparently. We don't do any modifications here, but now because this process is BuyNewShirt process, happens after the GetNewJob process, GetNewJob process, modifies the weight property and now BuyNewShirt is changed. Now the t-shirt size is a medium. I'll just give people a second to look at that, just digest it.

[4:54](https://www.youtube.com/watch?v=A357C-ux6Dw&t=294)
The core message that I want to get across here is that modifying the object, modifies it across all paths in the tag, the directed acyclic graph, the graph of processes, which can lead to complicated and time dependent bugs. We were seeing some problems in nf-core pipelines that would only appear if some processes took a little bit longer. That's what I was trying to demonstrate with that Browse process. They become really complicated, time dependent bugs to track down. Okay, so that's what the problem is. What is the solution? The solution is you should always return or whenever possible, return a new object instead of modifying objects in place. I am just going to introduce two very handy methods in Groovy for doing this, the most common modifications we do in metamap objects in nf-core. But the the idea that you should always return new objects is the general solution for this problem.

[6:05](https://www.youtube.com/watch?v=A357C-ux6Dw&t=365)
This is what the process was before, the process that I showed you before that had the bug. Modifying the property in place, but we can do this instead. Instead of returning the me object or the meta object, that map object transparently, oh, there's a bug here, that should be me and me, what we do is we create a new object and we add these two together and return that new object. This will fix the bug. It's important to note here that this plus operator is an alias for the .plus() method on maps in Groovy and the .plus() method will return a new map with precedence being given to the map on the right. This is the correct way. So this .plus() method is really important, it's a way of merging maps together in groovy. This is what the Groovy doc looks like with the link to the bottom. It returns a new map containing all the entries from the left and right, giving precedence to the right. That giving precedence to the right is important because it allows us to overwrite properties by placing them on the map on the right-hand side of the plus operator. Critically, the return object is a new map containing all the entries from the left and right.

[7:37](https://www.youtube.com/watch?v=A357C-ux6Dw&t=457)
Similarly, the same problem happens inside of map blocks or map closures. This is actually the more common case in the nf-core pipelines. What we're doing here is we're making the same object, a map of two keys named weight. Here I'm overriding the weight property or the weight key in this map by adding a new map. This will return a new map and make it clean for downstream use. In addition to overwrite, using this plus operator to overwrite properties in this map, I can also use it to overwrite properties and add new properties. Here I'm adding a new key, the employer key, and as Maxime said, I works at Seqera Labs in addition to overriding a property. It's a really nice way of piecing together and adding and overwriting maps in Nextflow pipelines.

[8:37](https://www.youtube.com/watch?v=A357C-ux6Dw&t=517)
The inverses operation you might think of is subtracting keys. This is also something that needs to happen quite a lot in nf-core pipelines where you want to take a subset of the keys. Rather than defining a new map, you can use this submap method. The submap method takes a collection of keys and returns a new map just containing those keys. Here we have a very verbose map, first name, last name, location, age, and employer. Let's say for downstream processes, I really only need first name and location. I can use this submap method, which does return a new object, doesn't modify the object in place, it returns a new object. You can see here in the return value in the documentation in the Groovy doc. That will be safe for modification. One complicated example in the wild here was from the Sarek pipeline. This was a really tricky bug to spot. It is now absolutely fixed in Sarek and is also fixed on the Nextflow side. But here we were taking the output of the FASTP process and we're taking the meta and the reads. We are calling the .sort() method on the reads. This is a dangerous operation because .sort() actually modifies the object in place. It sorts the object in place. Even though we're assigning it to a new value here, .sort() actually modifies reads, which had complicated implications for the publishing of those files. That's what you want to avoid. But for the most part, 99% of cases can be avoided by simply using the .plus() method on map objects or the .submap() method on objects for expanding and contracting meta map objects in nf-core pipelines.

[10:30](https://www.youtube.com/watch?v=A357C-ux6Dw&t=630)
Today's talk is very simple. Just has one clear message. Never modify anything in place and instead always return new objects when passing objects through nf-core and Nextflow pipelines. Have there any questions about that? I also have a VS code and we can go through examples interactively if people have more questions.

[10:55](https://www.youtube.com/watch?v=A357C-ux6Dw&t=655)
(host) Thank you very much, Rob, that was amazing and yes, I'm hoping people, does anyone have any questions here?

(question) I can always come up with a question. This is super insightful and I am sure I have written a lot of code which falls into these traps. How can we spot them? How can we spot these bugs?

(answer) It's almost always... let me just share my screen again. It's almost always this dot notation for modifying. But you can force it to be a little bit more clear. Actually, no, I'm not going to share my screen because I don't have an example of that. But you can force it. We could do one of two things in nf-core pipelines. We could instead of passing a map object, we could develop our own class that inherits from map and simply make it immutable. That would force things to, it wouldn't compile or not compile, it wouldn't run when you try to assign objects, try to modify it in place. That's one object. That's one path. That's a more secure path, but it does make things a little less convenient because there are sometimes, yeah, when you're creating an object where it's nice. But that might be one, and I'm toying around with what it would take to do that. Because I think it might be possible just to define a class in the lib directory, like a new metamap.groovy, change one invocation at the start of the pipeline and leave everything else unchanged. That might be the way we end up going, but I'm willing to give it a go and just see what people think.

(question cont.) I don't know how aware people are of the project, but I've just come off a call with Mattias and Julia, and we were talking about the nf-validation plugin, which Julia and Nicholas are writing, which has got a sample sheet, a new channel factory for generating a channel output from a sample sheet. That's where the metamap is most likely going to be created, right? At the point where we're passing the sample sheet. It could be something we could do within that Nextflow plugin. Then all the stuff, we don't have to touch the lib directory within the pipeline or anything. We could generate and define that class within the Nextflow plugin and hide it out of sight.

(answer cont.) Yeah. Yep. That's perfect. Yeah. That would work. Yeah. It's great.

[13:46](https://www.youtube.com/watch?v=A357C-ux6Dw&t=826)
(question) TZ asked a question about what was the best way to drop just one key from the metamap. My suggestion was to use minus and then submit with just that one key.

(answer) Yep. That's one option. Or you can pass in... Another option would be a .submap() and then it.keys() minus the key you don't want. That .submap() will take a collection of any kind. If you take a map and call the .keys() method, that returns a set of all the keys. You could subtract the key you don't want from that set and pass that it.keys() minus the key you don't want. But also the way that you've described there, Maxime, is also perfectly fine. I think there'd be about the same number of characters anyway. Personal preference.

(host) Okay. Thank you.

(answer continued) Yeah. There isn't a specific method for removing a key. But if we did create our own class, we could make our own method to remove a key. That might be enough.

(host) Yeah. But I think creating our own class would be something for another bytesize.

(speaker) Yeah. Yeah. Definitely.

[14:57](https://www.youtube.com/watch?v=A357C-ux6Dw&t=897)
(question) Is there a good source of documentation to read up on this anyway?

(answer) This idea is a little Nextflow based, so... no. But the Groovy docs that I linked in the slide. I'll pass a link on the Slack at the best place for Groovy docs.

(question cont.) I wonder if there's, I've forgotten who it was now, is it the Midnighter doc? Someone in the community has built a website, like a mate.materials website with common pitfalls and things for Nextflow. I wonder if it'd be a good one to go into that site.

(answer cont.) Oh yeah. I can suggest to add it there.

(question cont.) Apologies for that. I'm sure I'll be corrected in a second.

(comment) Yeah. I think it was Moritz. Moritz Viva. But also we've got some of this stuff in the advanced training docs, right?

(speaker) Yeah. Yeah. This is like a module in the advanced training. Yeah. We'll end up in training.sequera or training.nextflow.

[15:58](https://www.youtube.com/watch?v=A357C-ux6Dw&t=958)
(question) I just had a question. If you did want to sort, say reads as part of an operation for whatever reason, without actually appending or taking away from the map, what's the safest way then to just create a new copy with news.flow and would you reassign it? I mean, is reassigning safe? Is it just a pointer in memory then where you would be updating the original map anyway, you can just say new map is equal to old map, for example.

(answer) You can use... it requires... it's a really simple fix. All it requires is in .sort(), you pass it a TRUE as the first argument and that will sort it instead of sorting in place, but return a new object and then you can reassign it.

(question cont.) But is there a pointer in memory?

(answer cont.) I'm not sure I know what you mean.

(question cont.) If you create a new map, but this can be an issue in some programming languages where you have an individual map, right? So old map, and say you want to create a copy of it, you'd say new map is equal to old map, right? And then you do all of your downstream operations on new map. That's still change old map in place because it's just a pointer in memory to the old map. See what I mean?

(answer cont.) Yeah. If you pass .sort() TRUE, it will return a new map. The new map is a new object. The elements of the map will still be pointed to the same original copies, but that's okay because the order that you want to change. You'll get a new map and that object will have a new address in memory.

[17:45](https://www.youtube.com/watch?v=A357C-ux6Dw&t=1065)
(host) Okay. I think we have some more time. Does anyone have any more questions?

(speaker) Oh, Rike has just said she didn't know about exact(). That exact() method that I used as a convenient methods. The important thing to know about exact() is that it happens on the Nextflow head node rather than being farmed out to a process, which can be useful because particularly if you're operating in the cloud, you don't have to wait for VM speedups. It's also good for demos like this, because you can just write arbitrary Groovy. It's not that helpful. It's like in Nextflow pipelines, more often than not you don't need it. Sorry, Phil.

(comment) I was going to say, I usually am used to going as far as saying, don't use it because it's quite easy to abuse it and then crash the main Nextflow job or take down the head node.

(speaker) Yeah. Good for demos and maybe not best practice.

(host) Okay. Then I think we are good, so I will stop the recording.

</details>
