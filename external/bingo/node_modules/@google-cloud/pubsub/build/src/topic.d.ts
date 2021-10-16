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
/// <reference types="node" />
import { CallOptions } from 'google-gax';
import { google } from '../proto/pubsub';
import { IAM } from './iam';
import { Attributes, PublishCallback, Publisher, PublishOptions, PubsubMessage } from './publisher';
import { EmptyCallback, EmptyResponse, ExistsCallback, ExistsResponse, ObjectStream, PagedResponse, PageOptions, PubSub, RequestCallback, ResourceCallback } from './pubsub';
import { CreateSubscriptionCallback, CreateSubscriptionOptions, CreateSubscriptionResponse, Subscription, SubscriptionOptions } from './subscription';
export declare type TopicMetadata = google.pubsub.v1.ITopic;
declare type TopicCallback = ResourceCallback<Topic, TopicMetadata>;
declare type TopicResponse = [Topic, TopicMetadata];
export declare type CreateTopicCallback = TopicCallback;
export declare type CreateTopicResponse = TopicResponse;
export declare type GetTopicCallback = TopicCallback;
export declare type GetTopicResponse = TopicResponse;
export declare type GetTopicOptions = CallOptions & {
    autoCreate?: boolean;
};
declare type MetadataCallback = RequestCallback<TopicMetadata>;
declare type MetadataResponse = [TopicMetadata];
export declare type GetTopicMetadataCallback = MetadataCallback;
export declare type GetTopicMetadataResponse = MetadataResponse;
export declare type SetTopicMetadataCallback = MetadataCallback;
export declare type SetTopicMetadataResponse = MetadataResponse;
export declare type GetTopicSubscriptionsCallback = RequestCallback<Subscription, google.pubsub.v1.IListTopicSubscriptionsResponse>;
export declare type GetTopicSubscriptionsResponse = PagedResponse<Subscription, google.pubsub.v1.IListTopicSubscriptionsResponse>;
export declare type MessageOptions = PubsubMessage & {
    json?: any;
};
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
export declare class Topic {
    Promise?: PromiseConstructor;
    name: string;
    parent: PubSub;
    pubsub: PubSub;
    request: typeof PubSub.prototype.request;
    iam: IAM;
    metadata?: TopicMetadata;
    publisher: Publisher;
    getSubscriptionsStream: () => ObjectStream<Subscription>;
    constructor(pubsub: PubSub, name: string, options?: PublishOptions);
    create(gaxOpts?: CallOptions): Promise<CreateTopicResponse>;
    create(callback: CreateTopicCallback): void;
    create(gaxOpts: CallOptions, callback: CreateTopicCallback): void;
    createSubscription(name: string, callback: CreateSubscriptionCallback): void;
    createSubscription(name: string, options?: CreateSubscriptionOptions): Promise<CreateSubscriptionResponse>;
    createSubscription(name: string, options: CreateSubscriptionOptions, callback: CreateSubscriptionCallback): void;
    delete(callback: EmptyCallback): void;
    delete(gaxOpts?: CallOptions): Promise<EmptyResponse>;
    delete(gaxOpts: CallOptions, callback: EmptyCallback): void;
    exists(): Promise<ExistsResponse>;
    exists(callback: ExistsCallback): void;
    get(callback: GetTopicCallback): void;
    get(gaxOpts?: GetTopicOptions): Promise<GetTopicResponse>;
    get(gaxOpts: GetTopicOptions, callback: GetTopicCallback): void;
    getMetadata(callback: GetTopicMetadataCallback): void;
    getMetadata(gaxOpts: CallOptions, callback: GetTopicMetadataCallback): void;
    getMetadata(gaxOpts?: CallOptions): Promise<GetTopicMetadataResponse>;
    getSubscriptions(callback: GetTopicSubscriptionsCallback): void;
    getSubscriptions(options: PageOptions, callback: GetTopicSubscriptionsCallback): void;
    getSubscriptions(options?: PageOptions): Promise<GetTopicSubscriptionsResponse>;
    publish(data: Buffer, attributes?: Attributes): Promise<string>;
    publish(data: Buffer, callback: PublishCallback): void;
    publish(data: Buffer, attributes: Attributes, callback: PublishCallback): void;
    publishJSON(json: object, attributes?: Attributes): Promise<string>;
    publishJSON(json: object, callback: PublishCallback): void;
    publishJSON(json: object, attributes: Attributes, callback: PublishCallback): void;
    publishMessage(message: MessageOptions): Promise<[string]>;
    publishMessage(message: MessageOptions, callback: PublishCallback): void;
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
    resumePublishing(orderingKey: string): void;
    setMetadata(options: TopicMetadata, gaxOpts?: CallOptions): Promise<SetTopicMetadataResponse>;
    setMetadata(options: TopicMetadata, callback: SetTopicMetadataCallback): void;
    setMetadata(options: TopicMetadata, gaxOpts: CallOptions, callback: SetTopicMetadataCallback): void;
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
    setPublishOptions(options: PublishOptions): void;
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
    subscription(name: string, options?: SubscriptionOptions): Subscription;
    /**
     * Format the name of a topic. A Topic's full name is in the format of
     * 'projects/{projectId}/topics/{topicName}'.
     *
     * @private
     *
     * @return {string}
     */
    static formatName_(projectId: string, name: string): string;
}
export { PublishOptions };
