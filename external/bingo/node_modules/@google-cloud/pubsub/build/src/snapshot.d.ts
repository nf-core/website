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
import { CallOptions } from 'google-gax';
import { google } from '../proto/pubsub';
import { PubSub } from './index';
import { EmptyCallback, EmptyResponse, RequestCallback, ResourceCallback } from './pubsub';
import { Subscription } from './subscription';
export declare type CreateSnapshotCallback = ResourceCallback<Snapshot, google.pubsub.v1.ISnapshot>;
export declare type CreateSnapshotResponse = [Snapshot, google.pubsub.v1.ISnapshot];
export declare type SeekCallback = RequestCallback<google.pubsub.v1.ISeekResponse>;
export declare type SeekResponse = [google.pubsub.v1.ISeekResponse];
/**
 * A Snapshot object will give you access to your Cloud Pub/Sub snapshot.
 *
 * Snapshots are sometimes retrieved when using various methods:
 *
 * - {@link PubSub#getSnapshots}
 * - {@link PubSub#getSnapshotsStream}
 * - {@link PubSub#snapshot}
 *
 * Snapshots may be created with:
 *
 * - {@link Subscription#createSnapshot}
 *
 * You can use snapshots to seek a subscription to a specific point in time.
 *
 * - {@link Subscription#seek}
 *
 * @class
 *
 * @example
 * //-
 * // From {@link PubSub#getSnapshots}:
 * //-
 * pubsub.getSnapshots((err, snapshots) => {
 *   // `snapshots` is an array of Snapshot objects.
 * });
 *
 * //-
 * // From {@link PubSub#getSnapshotsStream}:
 * //-
 * pubsub.getSnapshotsStream()
 *   .on('error', console.error)
 *   .on('data', (snapshot) => {
 *     // `snapshot` is a Snapshot object.
 *   });
 *
 * //-
 * // From {@link PubSub#snapshot}:
 * //-
 * const snapshot = pubsub.snapshot('my-snapshot');
 * // snapshot is a Snapshot object.
 *
 * //-
 * // Create a snapshot with {module:pubsub/subscription#createSnapshot}:
 * //-
 * const subscription = pubsub.subscription('my-subscription');
 *
 * subscription.createSnapshot('my-snapshot', (err, snapshot) => {
 *   if (!err) {
 *     // `snapshot` is a Snapshot object.
 *   }
 * });
 *
 * //-
 * // Seek to your snapshot:
 * //-
 * const subscription = pubsub.subscription('my-subscription');
 *
 * subscription.seek('my-snapshot', (err) => {
 *   if (err) {
 *     // Error handling omitted.
 *   }
 * });
 */
export declare class Snapshot {
    parent: Subscription | PubSub;
    name: string;
    Promise?: PromiseConstructor;
    metadata?: google.pubsub.v1.ISnapshot;
    constructor(parent: Subscription | PubSub, name: string);
    delete(): Promise<EmptyResponse>;
    delete(callback: EmptyCallback): void;
    static formatName_(projectId: string, name: string): string;
    create(gaxOpts?: CallOptions): Promise<CreateSnapshotResponse>;
    create(callback: CreateSnapshotCallback): void;
    create(gaxOpts: CallOptions, callback: CreateSnapshotCallback): void;
    seek(gaxOpts?: CallOptions): Promise<SeekResponse>;
    seek(callback: SeekCallback): void;
    seek(gaxOpts: CallOptions, callback: SeekCallback): void;
}
