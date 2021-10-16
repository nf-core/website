/*!
 * Copyright 2019 Google Inc. All Rights Reserved.
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
import { BatchPublishOptions } from './message-batch';
import { Queue, OrderedQueue } from './message-queues';
import { Topic } from '../topic';
import { RequestCallback } from '../pubsub';
import { google } from '../../proto/pubsub';
export declare type PubsubMessage = google.pubsub.v1.IPubsubMessage;
export interface Attributes {
    [key: string]: string;
}
export declare type PublishCallback = RequestCallback<string>;
export interface PublishOptions {
    batching?: BatchPublishOptions;
    gaxOpts?: CallOptions;
    messageOrdering?: boolean;
}
export declare const BATCH_LIMITS: BatchPublishOptions;
/**
 * A Publisher object allows you to publish messages to a specific topic.
 *
 * @private
 * @class
 *
 * @see [Topics: publish API Documentation]{@link https://cloud.google.com/pubsub/docs/reference/rest/v1/projects.topics/publish}
 *
 * @param {Topic} topic The topic associated with this publisher.
 * @param {PublishOptions} [options] Configuration object.
 */
export declare class Publisher {
    Promise?: PromiseConstructor;
    topic: Topic;
    settings: PublishOptions;
    queue: Queue;
    orderedQueues: Map<string, OrderedQueue>;
    constructor(topic: Topic, options?: PublishOptions);
    publish(data: Buffer, attributes?: Attributes): Promise<string>;
    publish(data: Buffer, callback: PublishCallback): void;
    publish(data: Buffer, attributes: Attributes, callback: PublishCallback): void;
    /**
     * Publish the provided message.
     *
     * @private
     *
     * @throws {TypeError} If data is not a Buffer object.
     * @throws {TypeError} If any value in `attributes` object is not a string.
     *
     * @param {PubsubMessage} [message] Options for this message.
     * @param {PublishCallback} [callback] Callback function.
     */
    publishMessage(message: PubsubMessage, callback: PublishCallback): void;
    /**
     * Indicates to the publisher that it is safe to continue publishing for the
     * supplied ordering key.
     *
     * @private
     *
     * @param {string} key The ordering key to continue publishing for.
     */
    resumePublishing(key: string): void;
    /**
     * Sets the Publisher options.
     *
     * @private
     *
     * @param {PublishOptions} options The publisher options.
     */
    setOptions(options?: PublishOptions): void;
}
