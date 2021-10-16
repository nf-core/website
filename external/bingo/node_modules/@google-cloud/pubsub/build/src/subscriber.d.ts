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
import { PreciseDate } from '@google-cloud/precise-date';
import { EventEmitter } from 'events';
import { ClientStub } from 'google-gax';
import { google } from '../proto/pubsub';
import { FlowControlOptions } from './lease-manager';
import { BatchOptions } from './message-queues';
import { MessageStreamOptions } from './message-stream';
import { Subscription } from './subscription';
export declare type PullResponse = google.pubsub.v1.IPullResponse;
/**
 * Date object with nanosecond precision. Supports all standard Date arguments
 * in addition to several custom types.
 *
 * @external PreciseDate
 * @see {@link https://github.com/googleapis/nodejs-precise-date|PreciseDate}
 */
/**
 * Message objects provide a simple interface for users to get message data and
 * acknowledge the message.
 *
 * @example
 * subscription.on('message', message => {
 *   // {
 *   //   ackId: 'RUFeQBJMJAxESVMrQwsqWBFOBCEhPjA',
 *   //   attributes: {key: 'value'},
 *   //   data: Buffer.from('Hello, world!),
 *   //   id: '1551297743043',
 *   //   orderingKey: 'ordering-key',
 *   //   publishTime: new PreciseDate('2019-02-27T20:02:19.029534186Z'),
 *   //   received: 1551297743043,
 *   //   length: 13
 *   // }
 * });
 */
export declare class Message {
    ackId: string;
    attributes: {
        [key: string]: string;
    };
    data: Buffer;
    deliveryAttempt: number;
    id: string;
    orderingKey?: string;
    publishTime: PreciseDate;
    received: number;
    private _handled;
    private _length;
    private _subscriber;
    /**
     * @hideconstructor
     *
     * @param {Subscriber} sub The parent subscriber.
     * @param {object} message The raw message response.
     */
    constructor(sub: Subscriber, { ackId, message, deliveryAttempt }: google.pubsub.v1.IReceivedMessage);
    /**
     * The length of the message data.
     *
     * @type {number}
     */
    readonly length: number;
    /**
     * Acknowledges the message.
     *
     * @example
     * subscription.on('message', message => {
     *   message.ack();
     * });
     */
    ack(): void;
    /**
     * Modifies the ack deadline.
     *
     * @param {number} deadline The number of seconds to extend the deadline.
     * @private
     */
    modAck(deadline: number): void;
    /**
     * Removes the message from our inventory and schedules it to be redelivered.
     *
     * @example
     * subscription.on('message', message => {
     *   message.nack();
     * });
     */
    nack(): void;
}
export interface SubscriberOptions {
    ackDeadline?: number;
    batching?: BatchOptions;
    flowControl?: FlowControlOptions;
    streamingOptions?: MessageStreamOptions;
}
/**
 * @typedef {object} SubscriberOptions
 * @property {number} [ackDeadline=10] Acknowledge deadline in seconds. If left
 *     unset the initial value will be 10 seconds, but it will evolve into the
 *     99th percentile time it takes to acknowledge a message.
 * @property {BatchOptions} [batching] Request batching options.
 * @property {FlowControlOptions} [flowControl] Flow control options.
 * @property {MessageStreamOptions} [streamingOptions] Streaming options.
 */
/**
 * Subscriber class is used to manage all message related functionality.
 *
 * @private
 * @class
 *
 * @param {Subscription} subscription The corresponding subscription.
 * @param {SubscriberOptions} options The subscriber options.
 */
export declare class Subscriber extends EventEmitter {
    ackDeadline: number;
    isOpen: boolean;
    private _acks;
    private _histogram;
    private _inventory;
    private _isUserSetDeadline;
    private _latencies;
    private _modAcks;
    private _name;
    private _options;
    private _stream;
    private _subscription;
    constructor(subscription: Subscription, options?: {});
    /**
     * The 99th percentile of request latencies.
     *
     * @type {number}
     * @private
     */
    readonly modAckLatency: number;
    /**
     * The full name of the Subscription.
     *
     * @type {string}
     * @private
     */
    readonly name: string;
    /**
     * Acknowledges the supplied message.
     *
     * @param {Message} message The message to acknowledge.
     * @returns {Promise}
     * @private
     */
    ack(message: Message): Promise<void>;
    /**
     * Closes the subscriber. The returned promise will resolve once any pending
     * acks/modAcks are finished.
     *
     * @returns {Promise}
     * @private
     */
    close(): Promise<void>;
    /**
     * Gets the subscriber client instance.
     *
     * @returns {Promise<object>}
     * @private
     */
    getClient(): Promise<ClientStub>;
    /**
     * Modifies the acknowledge deadline for the provided message.
     *
     * @param {Message} message The message to modify.
     * @param {number} deadline The deadline.
     * @returns {Promise}
     * @private
     */
    modAck(message: Message, deadline: number): Promise<void>;
    /**
     * Modfies the acknowledge deadline for the provided message and then removes
     * it from our inventory.
     *
     * @param {Message} message The message.
     * @return {Promise}
     * @private
     */
    nack(message: Message): Promise<void>;
    /**
     * Starts pulling messages.
     * @private
     */
    open(): void;
    /**
     * Sets subscriber options.
     *
     * @param {SubscriberOptions} options The options.
     * @private
     */
    setOptions(options: SubscriberOptions): void;
    /**
     * Callback to be invoked when a new message is available.
     *
     * New messages will be added to the subscribers inventory, which in turn will
     * automatically extend the messages ack deadline until either:
     *   a. the user acks/nacks it
     *   b. the maxExtension option is hit
     *
     * If the message puts us at/over capacity, then we'll pause our message
     * stream until we've freed up some inventory space.
     *
     * New messages must immediately issue a ModifyAckDeadline request
     * (aka receipt) to confirm with the backend that we did infact receive the
     * message and its ok to start ticking down on the deadline.
     *
     * @private
     */
    private _onData;
    /**
     * Returns a promise that will resolve once all pending requests have settled.
     *
     * @private
     *
     * @returns {Promise}
     */
    private _waitForFlush;
}
