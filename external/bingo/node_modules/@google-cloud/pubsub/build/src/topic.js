"use strict";
/*!
 * Copyright 2017 Google Inc. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
Object.defineProperty(exports, "__esModule", { value: true });
const paginator_1 = require("@google-cloud/paginator");
const promisify_1 = require("@google-cloud/promisify");
const iam_1 = require("./iam");
const publisher_1 = require("./publisher");
/**
 * A Topic object allows you to interact with a Cloud Pub/Sub topic.
 *
 * @class
 * @param {PubSub} pubsub PubSub object.
 * @param {string} name Name of the topic.
 * @param {PublishOptions} [options] Publisher configuration object.
 *
 * @example
 * const {PubSub} = require('@google-cloud/pubsub');
 * const pubsub = new PubSub();
 *
 * const topic = pubsub.topic('my-topic');
 *
 * @example <caption>To enable message ordering, set `enableMessageOrdering` to true. Please note that this does not persist to an actual topic.</caption>
 * const topic = pubsub.topic('ordered-topic', {enableMessageOrdering: true});
 */
class Topic {
    constructor(pubsub, name, options) {
        this.getSubscriptionsStream = paginator_1.paginator.streamify('getSubscriptions');
        if (pubsub.Promise) {
            this.Promise = pubsub.Promise;
        }
        /**
         * The fully qualified name of this topic.
         * @name Topic#name
         * @type {string}
         */
        this.name = Topic.formatName_(pubsub.projectId, name);
        this.publisher = new publisher_1.Publisher(this, options);
        /**
         * The parent {@link PubSub} instance of this topic instance.
         * @name Topic#pubsub
         * @type {PubSub}
         */
        /**
         * The parent {@link PubSub} instance of this topic instance.
         * @name Topic#parent
         * @type {PubSub}
         */
        this.parent = this.pubsub = pubsub;
        this.request = pubsub.request.bind(pubsub);
        /**
         * [IAM (Identity and Access
         * Management)](https://cloud.google.com/pubsub/access_control) allows you
         * to set permissions on individual resources and offers a wider range of
         * roles: editor, owner, publisher, subscriber, and viewer. This gives you
         * greater flexibility and allows you to set more fine-grained access
         * control.
         *
         * *The IAM access control features described in this document are Beta,
         * including the API methods to get and set IAM policies, and to test IAM
         * permissions. Cloud Pub/Sub's use of IAM features is not covered by
         * any SLA or deprecation policy, and may be subject to
         * backward-incompatible changes.*
         *
         * @name Topic#iam
         * @mixes IAM
         *
         * @see [Access Control Overview]{@link https://cloud.google.com/pubsub/access_control}
         * @see [What is Cloud IAM?]{@link https://cloud.google.com/iam/}
         *
         * @example
         * const {PubSub} = require('@google-cloud/pubsub');
         * const pubsub = new PubSub();
         *
         * const topic = pubsub.topic('my-topic');
         *
         * //-
         * // Get the IAM policy for your topic.
         * //-
         * topic.iam.getPolicy((err, policy) => {
         *   console.log(policy);
         * });
         *
         * //-
         * // If the callback is omitted, we'll return a Promise.
         * //-
         * topic.iam.getPolicy().then((data) => {
         *   const policy = data[0];
         *   const apiResponse = data[1];
         * });
         */
        this.iam = new iam_1.IAM(pubsub, this.name);
    }
    /**
     * Create a topic.
     *
     * @param {object} [gaxOpts] Request configuration options, outlined
     *     here: https://googleapis.github.io/gax-nodejs/CallSettings.html.
     * @param {CreateTopicCallback} [callback] Callback function.
     * @returns {Promise<CreateTopicResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.create((err, topic, apiResponse) => {
     *   if (!err) {
     *     // The topic was created successfully.
     *   }
     * });
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.create().then((data) => {
     *   const topic = data[0];
     *   const apiResponse = data[1];
     * });
     */
    create(optsOrCallback, callback) {
        const gaxOpts = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        this.pubsub.createTopic(this.name, gaxOpts, callback);
    }
    /**
     * Create a subscription to this topic.
     *
     * @see [Subscriptions: create API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.subscriptions/create}
     *
     * @throws {Error} If subscription name is omitted.
     *
     * @param {string} name The name of the subscription.
     * @param {CreateSubscriptionRequest} [options] See a
     *     [Subscription
     * resource](https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.subscriptions).
     * @param {CreateSubscriptionCallback} [callback] Callback function.
     * @returns {Promise<CreateSubscriptionResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     * const callback = function(err, subscription, apiResponse) {};
     *
     * // Without specifying any options.
     * topic.createSubscription('newMessages', callback);
     *
     * // With options.
     * topic.createSubscription('newMessages', {
     *   ackDeadlineSeconds: 90
     * }, callback);
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.createSubscription('newMessages').then((data) => {
     *   const subscription = data[0];
     *   const apiResponse = data[1];
     * });
     */
    createSubscription(name, optsOrCallback, callback) {
        const options = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        this.pubsub.createSubscription(this, name, options, callback);
    }
    /**
     * Delete the topic. This will not delete subscriptions to this topic.
     *
     * @see [Topics: delete API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.topics/delete}
     *
     * @param {object} [gaxOpts] Request configuration options, outlined
     *     here: https://googleapis.github.io/gax-nodejs/CallSettings.html.
     * @param {function} [callback] The callback function.
     * @param {?error} callback.err An error returned while making this
     *     request.
     * @param {object} callback.apiResponse Raw API response.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.delete((err, apiResponse) => {});
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.delete().then((data) => {
     *   const apiResponse = data[0];
     * });
     */
    delete(optsOrCallback, callback) {
        const gaxOpts = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        const reqOpts = {
            topic: this.name,
        };
        this.request({
            client: 'PublisherClient',
            method: 'deleteTopic',
            reqOpts,
            gaxOpts: gaxOpts,
        }, callback);
    }
    /**
     * @typedef {array} TopicExistsResponse
     * @property {boolean} 0 Whether the topic exists
     */
    /**
     * @callback TopicExistsCallback
     * @param {?Error} err Request error, if any.
     * @param {boolean} exists Whether the topic exists.
     */
    /**
     * Check if a topic exists.
     *
     * @param {TopicExistsCallback} [callback] Callback function.
     * @returns {Promise<TopicExistsResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.exists((err, exists) => {});
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.exists().then((data) => {
     *   const exists = data[0];
     * });
     */
    exists(callback) {
        this.getMetadata(err => {
            if (!err) {
                callback(null, true);
                return;
            }
            if (err.code === 5) {
                callback(null, false);
                return;
            }
            callback(err);
        });
    }
    /**
     * @typedef {array} GetTopicResponse
     * @property {Topic} 0 The {@link Topic}.
     * @property {object} 1 The full API response.
     */
    /**
     * @callback GetTopicCallback
     * @param {?Error} err Request error, if any.
     * @param {Topic} topic The {@link Topic}.
     * @param {object} apiResponse The full API response.
     */
    /**
     * Get a topic if it exists.
     *
     * @param {object} [gaxOpts] Request configuration options, outlined
     *     here: https://googleapis.github.io/gax-nodejs/CallSettings.html.
     * @param {boolean} [gaxOpts.autoCreate=false] Automatically create the topic
     *     does not already exist.
     * @param {GetTopicCallback} [callback] Callback function.
     * @returns {Promise<GetTopicResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.get((err, topic, apiResponse) => {
     *   // The `topic` data has been populated.
     * });
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.get().then((data) => {
     *   const topic = data[0];
     *   const apiResponse = data[1];
     * });
     */
    get(optsOrCallback, callback) {
        const gaxOpts = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        const autoCreate = !!gaxOpts.autoCreate;
        delete gaxOpts.autoCreate;
        this.getMetadata(gaxOpts, (err, apiResponse) => {
            if (!err) {
                callback(null, this, apiResponse);
                return;
            }
            if (err.code !== 5 || !autoCreate) {
                callback(err, null, apiResponse);
                return;
            }
            this.create(gaxOpts, callback);
        });
    }
    /**
     * @typedef {array} GetTopicMetadataResponse
     * @property {object} 0 The full API response.
     */
    /**
     * @callback GetTopicMetadataCallback
     * @param {?Error} err Request error, if any.
     * @param {object} apiResponse The full API response.
     */
    /**
     * Get the official representation of this topic from the API.
     *
     * @see [Topics: get API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.topics/get}
     *
     * @param {object} [gaxOpts] Request configuration options, outlined
     *     here: https://googleapis.github.io/gax-nodejs/CallSettings.html.
     * @param {GetTopicMetadataCallback} [callback] Callback function.
     * @returns {Promise<GetTopicMetadataResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.getMetadata((err, apiResponse) => {});
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.getMetadata().then((data) => {
     *   const apiResponse = data[0];
     * });
     */
    getMetadata(optsOrCallback, callback) {
        const gaxOpts = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        const reqOpts = {
            topic: this.name,
        };
        this.request({
            client: 'PublisherClient',
            method: 'getTopic',
            reqOpts,
            gaxOpts: gaxOpts,
        }, (err, apiResponse) => {
            if (!err) {
                this.metadata = apiResponse;
            }
            callback(err, apiResponse);
        });
    }
    /**
     * Get a list of the subscriptions registered to this topic. You may
     * optionally provide a query object as the first argument to customize the
     * response.
     *
     * Your provided callback will be invoked with an error object if an API error
     * occurred or an array of {module:pubsub/subscription} objects.
     *
     * @see [Subscriptions: list API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.topics.subscriptions/list}
     *
     * @param {GetSubscriptionsRequest} [query] Query object for listing subscriptions.
     * @param {GetSubscriptionsCallback} [callback] Callback function.
     * @returns {Promise<GetSubscriptionsResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.getSubscriptions((err, subscriptions) => {
     *   // subscriptions is an array of `Subscription` objects.
     * });
     *
     * // Customize the query.
     * topic.getSubscriptions({
     *   pageSize: 3
     * }, callback);
     *
     * //-
     * // If the callback is omitted, we'll return a Promise.
     * //-
     * topic.getSubscriptions().then((data) => {
     *   const subscriptions = data[0];
     * });
     */
    getSubscriptions(optsOrCallback, callback) {
        const options = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        const reqOpts = Object.assign({
            topic: this.name,
        }, options);
        delete reqOpts.gaxOpts;
        delete reqOpts.autoPaginate;
        const gaxOpts = Object.assign({
            autoPaginate: options.autoPaginate,
        }, options.gaxOpts);
        this.request({
            client: 'PublisherClient',
            method: 'listTopicSubscriptions',
            reqOpts,
            gaxOpts,
        }, (err, subNames, ...args) => {
            let subscriptions;
            if (subNames) {
                subscriptions = subNames.map(sub => this.subscription(sub));
            }
            callback(err, subscriptions, ...args);
        });
    }
    /**
     * Publish the provided message.
     *
     * @deprecated Please use {@link Topic#publishMessage}.
     *
     * @throws {TypeError} If data is not a Buffer object.
     * @throws {TypeError} If any value in `attributes` object is not a string.
     *
     * @param {buffer} data The message data. This must come in the form of a
     *     Buffer object.
     * @param {object.<string, string>} [attributes] Attributes for this message.
     * @param {PublishCallback} [callback] Callback function.
     * @returns {Promise<PublishResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     * const data = Buffer.from('Hello, world!');
     *
     * const callback = (err, messageId) => {
     *   if (err) {
     *     // Error handling omitted.
     *   }
     * };
     *
     * topic.publish(data, callback);
     *
     * @example <caption>Optionally you can provide an object containing attributes for the message. Note that all values in the object must be strings.</caption>
     * const attributes = {
     *   key: 'value'
     * };
     *
     * topic.publish(data, attributes, callback);
     *
     * @example <caption>If the callback is omitted, we'll return a Promise.</caption>
     * topic.publish(data).then((messageId) => {});
     */
    publish(data, attrsOrCb, callback) {
        const attributes = typeof attrsOrCb === 'object' ? attrsOrCb : {};
        callback = typeof attrsOrCb === 'function' ? attrsOrCb : callback;
        return this.publishMessage({ data, attributes }, callback);
    }
    /**
     * Publish the provided JSON. It should be noted that all messages published
     * are done so in the form of a Buffer. This is simply a convenience method
     * that will transform JSON into a Buffer before publishing.
     * {@link Subscription} objects will always return message data in the form of
     * a Buffer, so any JSON published will require manual deserialization.
     *
     * @deprecated Please use the `json` option via {@link Topic#publishMessage}.
     *
     * @throws {Error} If non-object data is provided.
     *
     * @param {object} json The JSON data to publish.
     * @param {object} [attributes] Attributes for this message.
     * @param {PublishCallback} [callback] Callback function.
     * @returns {Promise<PublishResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     * const topic = pubsub.topic('my-topic');
     *
     * const data = {
     *   foo: 'bar'
     * };
     *
     * const callback = (err, messageId) => {
     *   if (err) {
     *     // Error handling omitted.
     *   }
     * };
     *
     * topic.publishJSON(data, callback);
     *
     * @example <caption>Optionally you can provide an object containing attributes for the message. Note that all values in the object must be strings.</caption>
     * const attributes = {
     *   key: 'value'
     * };
     *
     * topic.publishJSON(data, attributes, callback);
     *
     * @example <caption>If the callback is omitted, we'll return a Promise.</caption>
     * topic.publishJSON(data).then((messageId) => {});
     */
    publishJSON(json, attrsOrCb, callback) {
        if (!json || typeof json !== 'object') {
            throw new Error('First parameter should be an object.');
        }
        const attributes = typeof attrsOrCb === 'object' ? attrsOrCb : {};
        callback = typeof attrsOrCb === 'function' ? attrsOrCb : callback;
        return this.publishMessage({ json, attributes }, callback);
    }
    /**
     * @typedef {object} MessageOptions
     * @property {buffer} [data] The message data.
     * @property {object} [json] Convenience property to publish JSON data. This
     *     will transform the provided JSON into a Buffer before publishing.
     *     {@link Subscription} objects will always return message data in the
     *     form of a Buffer, so any JSON published will require manual
     *     deserialization.
     * @property {object.<string, string>} [attributes] Attributes for this
     *     message.
     * @property {string} [orderingKey] A message ordering key.
     */
    /**
     * Publish the provided message.
     *
     * @throws {TypeError} If data is not a Buffer object.
     * @throws {TypeError} If any value in `attributes` object is not a string.
     *
     * @param {MessageOptions} message Message object.
     * @param {PublishCallback} [callback] Callback function.
     * @returns {Promise<PublishResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     * const topic = pubsub.topic('my-topic');
     *
     * const data = Buffer.from('Hello, world!');
     *
     * const callback = (err, messageId) => {
     *   if (err) {
     *     // Error handling omitted.
     *   }
     * };
     *
     * topic.publishMessage({data}, callback);
     *
     * @example <caption>Publish JSON message data.</caption>
     * const json = {foo: 'bar'};
     *
     * topic.publishMessage({json}, callback);
     *
     * @example <caption>To publish messages in order (this is still experimental), make sure message ordering is enabled and provide an ordering key</caption>
     * const topic = pubsub.topic('ordered-topic', {messageOrdering: true});
     * const orderingKey = 'my-key';
     *
     * topic.publishMessage({data, orderingKey}, callback);
     *
     * @example <caption>If the callback is omitted, we'll return a Promise.</caption>
     * const [messageId] = await topic.publishMessage({data});
     */
    publishMessage(message, callback) {
        message = Object.assign({}, message);
        if (message.json && typeof message.json === 'object') {
            message.data = Buffer.from(JSON.stringify(message.json));
            delete message.json;
        }
        return this.publisher.publishMessage(message, callback);
    }
    /**
     * In the event that the client fails to publish an ordered message, all
     * subsequent publish calls using the same ordering key will fail. Calling
     * this method will disregard the publish failure, allowing the supplied
     * ordering key to be used again in the future.
     *
     * @param {string} orderingKey The ordering key in question.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     * const topic = pubsub.topic('my-topic', {messageOrdering: true});
     *
     * const orderingKey = 'foo';
     * const data = Buffer.from('Hello, order!');
     *
     * topic.publishMessage({data, orderingKey}, err => {
     *   if (err) {
     *     topic.resumePublishing(orderingKey);
     *   }
     * });
     */
    resumePublishing(orderingKey) {
        this.publisher.resumePublishing(orderingKey);
    }
    /**
     * @typedef {array} SetTopicMetadataResponse
     * @property {object} 0 The full API response.
     */
    /**
     * @callback SetTopicMetadataCallback
     * @param {?Error} err Request error, if any.
     * @param {object} apiResponse The full API response.
     */
    /**
     * Updates the topic.
     *
     * @see [UpdateTopicRequest API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/UpdateTopicRequest}
     *
     * @param {object} metadata The fields to update. This should be structured
     *     like a {@link
     * https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.topics#Topic|Topic
     * object}.
     * @param {object} [gaxOpts] Request configuration options, outlined
     *     here: https://googleapis.github.io/gax-nodejs/CallSettings.html.
     * @param {SetTopicMetadataCallback} [callback] Callback function.
     * @returns {Promise<SetTopicMetadataResponse>}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     * const metadata = {
     *   labels: {foo: 'bar'}
     * };
     *
     * topic.setMetadata(metadata, err => {
     *   if (err) {
     *     // Error handling omitted.
     *   }
     * });
     *
     * @example <caption>If the callback is omitted, we'll return a
     * Promise.</caption>
     * topic.setMetadata(metadata).then((data) => {
     *   const apiResponse = data[0];
     * });
     */
    setMetadata(options, optsOrCallback, callback) {
        const gaxOpts = typeof optsOrCallback === 'object' ? optsOrCallback : {};
        callback = typeof optsOrCallback === 'function' ? optsOrCallback : callback;
        const topic = Object.assign({ name: this.name }, options);
        const updateMask = { paths: Object.keys(options) };
        const reqOpts = { topic, updateMask };
        this.request({
            client: 'PublisherClient',
            method: 'updateTopic',
            reqOpts,
            gaxOpts,
        }, callback);
    }
    /**
     * Set the publisher options.
     *
     * @param {PublishOptions} options The publisher options.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     *
     * topic.setPublishOptions({
     *   batching: {
     *     maxMilliseconds: 10
     *   }
     * });
     */
    setPublishOptions(options) {
        this.publisher.setOptions(options);
    }
    /**
     * Create a Subscription object. This command by itself will not run any API
     * requests. You will receive a {module:pubsub/subscription} object,
     * which will allow you to interact with a subscription.
     *
     * @throws {Error} If subscription name is omitted.
     *
     * @param {string} name Name of the subscription.
     * @param {SubscriberOptions} [options] Configuration object.
     * @return {Subscription}
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     * const subscription = topic.subscription('my-subscription');
     *
     * // Register a listener for `message` events.
     * subscription.on('message', (message) => {
     *   // Called every time a message is received.
     *   // message.id = ID of the message.
     *   // message.ackId = ID used to acknowledge the message receival.
     *   // message.data = Contents of the message.
     *   // message.attributes = Attributes of the message.
     *   // message.publishTime = Timestamp when Pub/Sub received the message.
     * });
     */
    subscription(name, options) {
        options = options || {};
        options.topic = this;
        return this.pubsub.subscription(name, options);
    }
    /**
     * Format the name of a topic. A Topic's full name is in the format of
     * 'projects/{projectId}/topics/{topicName}'.
     *
     * @private
     *
     * @return {string}
     */
    static formatName_(projectId, name) {
        // Simple check if the name is already formatted.
        if (name.indexOf('/') > -1) {
            return name;
        }
        return 'projects/' + projectId + '/topics/' + name;
    }
}
exports.Topic = Topic;
/**
 * Get a list of the {module:pubsub/subscription} objects registered to this
 * topic as a readable object stream.
 *
 * @method PubSub#getSubscriptionsStream
 * @param {GetSubscriptionsRequest} [options] Configuration object. See
 *     {@link PubSub#getSubscriptions} for a complete list of options.
 * @returns {ReadableStream} A readable stream of {@link Subscription} instances.
 *
 * @example
 * const {PubSub} = require('@google-cloud/pubsub');
 * const pubsub = new PubSub();
 *
 * const topic = pubsub.topic('my-topic');
 *
 * topic.getSubscriptionsStream()
 *   .on('error', console.error)
 *   .on('data', (subscription) => {
 *     // subscription is a Subscription object.
 *   })
 *   .on('end', () => {
 *     // All subscriptions retrieved.
 *   });
 *
 * //-
 * // If you anticipate many results, you can end a stream early to prevent
 * // unnecessary processing and API requests.
 * //-
 * topic.getSubscriptionsStream()
 *   .on('data', function(subscription) {
 *     this.end();
 *   });
 */
/*! Developer Documentation
 *
 * These methods can be agto-paginated.
 */
paginator_1.paginator.extend(Topic, ['getSubscriptions']);
/*! Developer Documentation
 *
 * All async methods (except for streams) will return a Promise in the event
 * that a callback is omitted.
 */
promisify_1.promisifyAll(Topic, {
    exclude: [
        'publish',
        'publishJSON',
        'publishMessage',
        'setPublishOptions',
        'subscription',
    ],
});
//# sourceMappingURL=topic.js.map