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
import { Metadata, ServiceError, status, StatusObject } from '@grpc/grpc-js';
import { PassThrough } from 'stream';
import { Subscriber } from './subscriber';
/**
 * Error wrapper for gRPC status objects.
 *
 * @class
 *
 * @param {object} status The gRPC status object.
 */
export declare class StatusError extends Error implements ServiceError {
    code: status;
    details: string;
    metadata: Metadata;
    constructor(status: StatusObject);
}
/**
 * Error thrown when we fail to open a channel for the message stream.
 *
 * @class
 *
 * @param {Error} err The original error.
 */
export declare class ChannelError extends Error implements ServiceError {
    code: status;
    details: string;
    metadata: Metadata;
    constructor(err: Error);
}
export interface MessageStreamOptions {
    highWaterMark?: number;
    maxStreams?: number;
    timeout?: number;
}
/**
 * @typedef {object} MessageStreamOptions
 * @property {number} [highWaterMark=0] Configures the Buffer level for all
 *     underlying streams. See
 *     {@link https://nodejs.org/en/docs/guides/backpressuring-in-streams/} for
 *     more details.
 * @property {number} [maxStreams=5] Number of streaming connections to make.
 * @property {number} [timeout=300000] Timeout for establishing a connection.
 */
/**
 * Streaming class used to manage multiple StreamingPull requests.
 *
 * @private
 * @class
 *
 * @param {Subscriber} sub The parent subscriber.
 * @param {MessageStreamOptions} [options] The message stream options.
 */
export declare class MessageStream extends PassThrough {
    destroyed: boolean;
    private _keepAliveHandle;
    private _fillHandle?;
    private _options;
    private _retrier;
    private _streams;
    private _subscriber;
    constructor(sub: Subscriber, options?: MessageStreamOptions);
    /**
     * Destroys the stream and any underlying streams.
     *
     * @param {error?} err An error to emit, if any.
     * @private
     */
    destroy(err?: Error): void;
    /**
     * Adds a StreamingPull stream to the combined stream.
     *
     * @private
     *
     * @param {stream} stream The StreamingPull stream.
     */
    private _addStream;
    /**
     * Attempts to create and cache the desired number of StreamingPull requests.
     * gRPC does not supply a way to confirm that a stream is connected, so our
     * best bet is to open the streams and use the client.waitForReady() method to
     * confirm everything is ok.
     *
     * @private
     *
     * @returns {Promise}
     */
    private _fillStreamPool;
    /**
     * It is critical that we keep as few `PullResponse` objects in memory as
     * possible to reduce the number of potential redeliveries. Because of this we
     * want to bypass gax for StreamingPull requests to avoid creating a Duplexify
     * stream, doing so essentially doubles the size of our readable buffer.
     *
     * @private
     *
     * @returns {Promise.<object>}
     */
    private _getClient;
    /**
     * Since we do not use the streams to ack/modAck messages, they will close
     * by themselves unless we periodically send empty messages.
     *
     * @private
     */
    private _keepAlive;
    /**
     * Once the stream has nothing left to read, we'll remove it and attempt to
     * refill our stream pool if needed.
     *
     * @private
     *
     * @param {Duplex} stream The ended stream.
     * @param {object} status The stream status.
     */
    private _onEnd;
    /**
     * gRPC will usually emit a status as a ServiceError via `error` event before
     * it emits the status itself. In order to cut back on emitted errors, we'll
     * wait a tick on error and ignore it if the status has been received.
     *
     * @private
     *
     * @param {stream} stream The stream that errored.
     * @param {Error} err The error.
     */
    private _onError;
    /**
     * gRPC streams will emit a status event once the connection has been
     * terminated. This is preferable to end/close events because we'll receive
     * information as to why the stream closed and if it is safe to open another.
     *
     * @private
     *
     * @param {stream} stream The stream that was closed.
     * @param {object} status The status message stating why it was closed.
     */
    private _onStatus;
    /**
     * Removes a stream from the combined stream.
     *
     * @private
     *
     * @param {stream} stream The stream to remove.
     */
    private _removeStream;
    /**
     * Neither gRPC or gax allow for the highWaterMark option to be specified.
     * However using the default value (16) it is possible to end up with a lot of
     * PullResponse objects stored in internal buffers. If this were to happen
     * and the client were slow to process messages, we could potentially see a
     * very large number of redeliveries happen before the messages even made it
     * to the client.
     *
     * @private
     *
     * @param {Duplex} stream The duplex stream to adjust the
     *     highWaterMarks for.
     */
    private _setHighWaterMark;
    /**
     * Promisified version of gRPCs Client#waitForReady function.
     *
     * @private
     *
     * @param {object} client The gRPC client to wait for.
     * @returns {Promise}
     */
    private _waitForClientReady;
}
