/*!
 * Copyright 2014 Google Inc. All Rights Reserved.
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
import { GoogleAuth } from 'google-auth-library';
import * as gax from 'google-gax';
import { ServiceError, ChannelCredentials } from '@grpc/grpc-js';
import { Snapshot } from './snapshot';
import { Subscription, SubscriptionOptions, CreateSubscriptionOptions, CreateSubscriptionCallback, CreateSubscriptionResponse } from './subscription';
import { Topic, GetTopicSubscriptionsCallback, GetTopicSubscriptionsResponse, CreateTopicCallback, CreateTopicResponse } from './topic';
import { PublishOptions } from './publisher';
import { CallOptions } from 'google-gax';
import { Transform } from 'stream';
import { google } from '../proto/pubsub';
export declare type Omit<T, K extends keyof T> = Pick<T, Exclude<keyof T, K>>;
export interface ClientConfig extends gax.GrpcClientOptions {
    apiEndpoint?: string;
    servicePath?: string;
    port?: string | number;
    sslCreds?: ChannelCredentials;
}
export interface PageOptions {
    gaxOpts?: CallOptions;
    pageSize?: number;
    pageToken?: string;
    autoPaginate?: boolean;
}
export declare type GetSnapshotsCallback = RequestCallback<Snapshot, google.pubsub.v1.IListSnapshotsResponse>;
export declare type GetSnapshotsResponse = PagedResponse<Snapshot, google.pubsub.v1.IListSnapshotsResponse>;
export declare type GetSubscriptionsOptions = PageOptions & {
    topic?: string | Topic;
};
declare type GetAllSubscriptionsCallback = RequestCallback<Subscription, google.pubsub.v1.IListSubscriptionsResponse>;
declare type GetAllSubscriptionsResponse = PagedResponse<Subscription, google.pubsub.v1.IListSubscriptionsResponse>;
export declare type GetSubscriptionsCallback = GetAllSubscriptionsCallback | GetTopicSubscriptionsCallback;
export declare type GetSubscriptionsResponse = GetAllSubscriptionsResponse | GetTopicSubscriptionsResponse;
export declare type GetTopicsCallback = RequestCallback<Topic, google.pubsub.v1.IListTopicsResponse>;
export declare type GetTopicsResponse = PagedResponse<Topic, google.pubsub.v1.IListTopicsResponse>;
export declare type EmptyCallback = RequestCallback<google.protobuf.IEmpty>;
export declare type EmptyResponse = [google.protobuf.IEmpty];
export declare type ExistsCallback = RequestCallback<boolean>;
export declare type ExistsResponse = [boolean];
export interface GetClientConfig {
    client: 'PublisherClient' | 'SubscriberClient';
}
export interface RequestConfig extends GetClientConfig {
    method: string;
    reqOpts?: object;
    gaxOpts?: CallOptions;
}
export interface ResourceCallback<Resource, Response> {
    (err: ServiceError | null, resource?: Resource | null, response?: Response | null): void;
}
export declare type RequestCallback<T, R = void> = R extends void ? NormalCallback<T> : PagedCallback<T, R>;
export interface NormalCallback<TResponse> {
    (err: ServiceError | null, res?: TResponse | null): void;
}
export interface PagedCallback<Item, Response> {
    (err: ServiceError | null, results?: Item[] | null, nextQuery?: {} | null, response?: Response | null): void;
}
export declare type PagedResponse<Item, Response> = [Item[]] | [Item[], {} | null, Response];
export declare type ObjectStream<O> = {
    addListener(event: 'data', listener: (data: O) => void): ObjectStream<O>;
    emit(event: 'data', data: O): boolean;
    on(event: 'data', listener: (data: O) => void): ObjectStream<O>;
    once(event: 'data', listener: (data: O) => void): ObjectStream<O>;
    prependListener(event: 'data', listener: (data: O) => void): ObjectStream<O>;
    prependOnceListener(event: 'data', listener: (data: O) => void): ObjectStream<O>;
} & Transform;
interface GetClientCallback {
    (err: Error | null, gaxClient?: gax.ClientStub): void;
}
/**
 * @typedef {object} ClientConfig
 * @property {string} [projectId] The project ID from the Google Developer's
 *     Console, e.g. 'grape-spaceship-123'. We will also check the environment
 *     variable `GCLOUD_PROJECT` for your project ID. If your app is running in
 *     an environment which supports {@link
 * https://cloud.google.com/docs/authentication/production#providing_credentials_to_your_application
 * Application Default Credentials}, your project ID will be detected
 * automatically.
 * @property {string} [keyFilename] Full path to the a .json, .pem, or .p12 key
 *     downloaded from the Google Developers Console. If you provide a path to a
 *     JSON file, the `projectId` option above is not necessary. NOTE: .pem and
 *     .p12 require you to specify the `email` option as well.
 * @property {string} [apiEndpoint] The `apiEndpoint` from options will set the
 *     host. If not set, the `PUBSUB_EMULATOR_HOST` environment variable from
 *     the gcloud SDK is honored, otherwise the actual API endpoint will be
 *     used.
 * @property {string} [email] Account email address. Required when using a .pem
 *     or .p12 keyFilename.
 * @property {object} [credentials] Credentials object.
 * @property {string} [credentials.client_email]
 * @property {string} [credentials.private_key]
 * @property {boolean} [autoRetry=true] Automatically retry requests if the
 *     response is related to rate limits or certain intermittent server errors.
 *     We will exponentially backoff subsequent requests by default.
 * @property {number} [maxRetries=3] Maximum number of automatic retries
 *     attempted before returning the error.
 * @property {Constructor} [promise] Custom promise module to use instead of
 *     native Promises.
 */
/**
 * [Cloud Pub/Sub](https://developers.google.com/pubsub/overview) is a
 * reliable, many-to-many, asynchronous messaging service from Cloud
 * Platform.
 *
 * @class
 *
 * @see [Cloud Pub/Sub overview]{@link https://developers.google.com/pubsub/overview}
 *
 * @param {ClientConfig} [options] Configuration options.
 *
 * @example <caption>Import the client library</caption>
 * const {PubSub} = require('@google-cloud/pubsub');
 *
 * @example <caption>Create a client that uses <a href="https://cloud.google.com/docs/authentication/production#providing_credentials_to_your_application">Application Default Credentials (ADC)</a>:</caption>
 * const pubsub = new PubSub();
 *
 * @example <caption>Create a client with <a href="https://cloud.google.com/docs/authentication/production#obtaining_and_providing_service_account_credentials_manually">explicit credentials</a>:</caption>
 * const pubsub = new PubSub({
 *   projectId: 'your-project-id',
 *   keyFilename: '/path/to/keyfile.json'
 * });
 *
 * @example <caption>include:samples/quickstart.js</caption>
 * region_tag:pubsub_quickstart_create_topic
 * Full quickstart example:
 */
export declare class PubSub {
    options: ClientConfig;
    isEmulator: boolean;
    api: {
        [key: string]: gax.ClientStub;
    };
    auth: GoogleAuth;
    projectId: string;
    Promise?: PromiseConstructor;
    getSubscriptionsStream: () => ObjectStream<Subscription>;
    getSnapshotsStream: () => ObjectStream<Snapshot>;
    getTopicsStream: () => ObjectStream<Topic>;
    constructor(options?: ClientConfig);
    createSubscription(topic: Topic | string, name: string, options?: CreateSubscriptionOptions): Promise<CreateSubscriptionResponse>;
    createSubscription(topic: Topic | string, name: string, callback: CreateSubscriptionCallback): void;
    createSubscription(topic: Topic | string, name: string, options: CreateSubscriptionOptions, callback: CreateSubscriptionCallback): void;
    createTopic(name: string, gaxOpts?: CallOptions): Promise<CreateTopicResponse>;
    createTopic(name: string, callback: CreateTopicCallback): void;
    createTopic(name: string, gaxOpts: CallOptions, callback: CreateTopicCallback): void;
    /**
     * Determine the appropriate endpoint to use for API requests, first trying
     * the local `apiEndpoint` parameter. If the `apiEndpoint` parameter is null
     * we try Pub/Sub emulator environment variable (PUBSUB_EMULATOR_HOST),
     * otherwise the default JSON API.
     *
     * @private
     */
    determineBaseUrl_(): void;
    getSnapshots(options?: PageOptions): Promise<GetSnapshotsResponse>;
    getSnapshots(callback: GetSnapshotsCallback): void;
    getSnapshots(options: PageOptions, callback: GetSnapshotsCallback): void;
    getSubscriptions(options?: GetSubscriptionsOptions): Promise<GetSubscriptionsResponse>;
    getSubscriptions(callback: GetSubscriptionsCallback): void;
    getSubscriptions(options: GetSubscriptionsOptions, callback: GetSubscriptionsCallback): void;
    getTopics(options?: PageOptions): Promise<GetTopicsResponse>;
    getTopics(callback: GetTopicsCallback): void;
    getTopics(options: PageOptions, callback: GetTopicsCallback): void;
    /**
     * Callback function to PubSub.getClient_().
     * @private
     * @callback GetClientCallback
     * @param err - Error, if any.
     * @param gaxClient - The gax client specified in RequestConfig.client.
     *                    Typed any since it's importing Javascript source.
     */
    /**
     * Get the PubSub client object.
     *
     * @private
     *
     * @param {object} config Configuration object.
     * @param {object} config.gaxOpts GAX options.
     * @param {function} config.method The gax method to call.
     * @param {object} config.reqOpts Request options.
     * @param {function} [callback] The callback function.
     */
    getClient_(config: GetClientConfig, callback: GetClientCallback): void;
    /**
     * Get the PubSub client object.
     *
     * @private
     *
     * @param {object} config Configuration object.
     * @param {object} config.gaxOpts GAX options.
     * @param {function} config.method The gax method to call.
     * @param {object} config.reqOpts Request options.
     * @returns {Promise}
     */
    getClientAsync_(config: GetClientConfig): Promise<gax.ClientStub>;
    /**
     * Funnel all API requests through this method, to be sure we have a project
     * ID.
     *
     * @private
     *
     * @param {object} config Configuration object.
     * @param {object} config.gaxOpts GAX options.
     * @param {function} config.method The gax method to call.
     * @param {object} config.reqOpts Request options.
     * @param {function} [callback] The callback function.
     */
    request<T, R = void>(config: RequestConfig, callback: RequestCallback<T, R>): void;
    /**
     * Create a Snapshot object. See {@link Subscription#createSnapshot} to
     * create a snapshot.
     *
     * @throws {Error} If a name is not provided.
     *
     * @param {string} name The name of the snapshot.
     * @returns {Snapshot} A {@link Snapshot} instance.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const snapshot = pubsub.snapshot('my-snapshot');
     */
    snapshot(name: string): Snapshot;
    /**
     * Create a Subscription object. This command by itself will not run any API
     * requests. You will receive a {@link Subscription} object,
     * which will allow you to interact with a subscription.
     *
     * @throws {Error} If subscription name is omitted.
     *
     * @param {string} name Name of the subscription.
     * @param {SubscriberOptions} [options] Configuration object.
     * @returns {Subscription} A {@link Subscription} instance.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const subscription = pubsub.subscription('my-subscription');
     *
     * // Register a listener for `message` events.
     * subscription.on('message', function(message) {
     *   // Called every time a message is received.
     *   // message.id = ID of the message.
     *   // message.ackId = ID used to acknowledge the message receival.
     *   // message.data = Contents of the message.
     *   // message.attributes = Attributes of the message.
     *   // message.publishTime = Date when Pub/Sub received the message.
     * });
     */
    subscription(name: string, options?: SubscriptionOptions): Subscription;
    /**
     * Create a Topic object. See {@link PubSub#createTopic} to create a topic.
     *
     * @throws {Error} If a name is not provided.
     *
     * @param {string} name The name of the topic.
     * @returns {Topic} A {@link Topic} instance.
     *
     * @example
     * const {PubSub} = require('@google-cloud/pubsub');
     * const pubsub = new PubSub();
     *
     * const topic = pubsub.topic('my-topic');
     */
    topic(name: string, options?: PublishOptions): Topic;
}
export {};
