/*!
 * Copyright 2018 Google Inc. All Rights Reserved.
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
import { Metadata, ServiceError, status } from '@grpc/grpc-js';
import defer = require('p-defer');
import { Message, Subscriber } from './subscriber';
declare type QueuedMessages = Array<[string, number?]>;
export interface BatchOptions {
    callOptions?: CallOptions;
    maxMessages?: number;
    maxMilliseconds?: number;
}
/**
 * Error class used to signal a batch failure.
 *
 * @class
 *
 * @param {string} message The error message.
 * @param {ServiceError} err The grpc service error.
 */
export declare class BatchError extends Error implements ServiceError {
    ackIds: string[];
    code: status;
    details: string;
    metadata: Metadata;
    constructor(err: ServiceError, ackIds: string[], rpc: string);
}
/**
 * @typedef {object} BatchOptions
 * @property {object} [callOptions] Request configuration option, outlined
 *     here: {@link https://googleapis.github.io/gax-nodejs/interfaces/CallOptions.html}.
 * @property {number} [maxMessages=3000] Maximum number of messages allowed in
 *     each batch sent.
 * @property {number} [maxMilliseconds=100] Maximum duration to wait before
 *     sending a batch. Batches can be sent earlier if the maxMessages option
 *     is met before the configured duration has passed.
 */
/**
 * Class for buffering ack/modAck requests.
 *
 * @private
 * @class
 *
 * @param {Subscriber} sub The subscriber we're queueing requests for.
 * @param {BatchOptions} options Batching options.
 */
export declare abstract class MessageQueue {
    numPendingRequests: number;
    protected _onFlush?: defer.DeferredPromise<void>;
    protected _options: BatchOptions;
    protected _requests: QueuedMessages;
    protected _subscriber: Subscriber;
    protected _timer?: NodeJS.Timer;
    protected abstract _sendBatch(batch: QueuedMessages): Promise<void>;
    constructor(sub: Subscriber, options?: BatchOptions);
    /**
     * Gets the default buffer time in ms.
     *
     * @returns {number}
     * @private
     */
    readonly maxMilliseconds: number;
    /**
     * Adds a message to the queue.
     *
     * @param {Message} message The message to add.
     * @param {number} [deadline] The deadline.
     * @private
     */
    add({ ackId }: Message, deadline?: number): void;
    /**
     * Sends a batch of messages.
     * @private
     */
    flush(): Promise<void>;
    /**
     * Returns a promise that resolves after the next flush occurs.
     *
     * @returns {Promise}
     * @private
     */
    onFlush(): Promise<void>;
    /**
     * Set the batching options.
     *
     * @param {BatchOptions} options Batching options.
     * @private
     */
    setOptions(options: BatchOptions): void;
}
/**
 * Queues up Acknowledge (ack) requests.
 *
 * @private
 * @class
 */
export declare class AckQueue extends MessageQueue {
    /**
     * Sends a batch of ack requests.
     *
     * @private
     *
     * @param {Array.<Array.<string|number>>} batch Array of ackIds and deadlines.
     * @return {Promise}
     */
    protected _sendBatch(batch: QueuedMessages): Promise<void>;
}
/**
 * Queues up ModifyAckDeadline requests and sends them out in batches.
 *
 * @private
 * @class
 */
export declare class ModAckQueue extends MessageQueue {
    /**
     * Sends a batch of modAck requests. Each deadline requires its own request,
     * so we have to group all the ackIds by deadline and send multiple requests.
     *
     * @private
     *
     * @param {Array.<Array.<string|number>>} batch Array of ackIds and deadlines.
     * @return {Promise}
     */
    protected _sendBatch(batch: QueuedMessages): Promise<void>;
}
export {};
