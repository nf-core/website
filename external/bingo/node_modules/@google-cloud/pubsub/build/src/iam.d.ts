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
import { CallOptions } from 'google-gax';
import { google } from '../proto/iam';
import { Omit, PubSub, RequestCallback, ResourceCallback } from './pubsub';
export declare type Policy = {
    etag?: string | Buffer;
} & Omit<google.iam.v1.IPolicy, 'etag'>;
export declare type GetPolicyCallback = RequestCallback<Policy>;
export declare type SetPolicyCallback = RequestCallback<Policy>;
export declare type SetPolicyResponse = [Policy];
export declare type GetPolicyResponse = [Policy];
/**
 * Shows which IAM permissions is allowed.
 * The key to this object are the IAM permissions (string) and the values are
 * booleans, true if permissions are granted to the corresponding key.
 */
export interface IamPermissionsMap {
    [key: string]: boolean;
}
export declare type TestIamPermissionsResponse = [IamPermissionsMap, google.iam.v1.ITestIamPermissionsResponse];
export declare type TestIamPermissionsCallback = ResourceCallback<IamPermissionsMap, google.iam.v1.ITestIamPermissionsResponse>;
/**
 * [IAM (Identity and Access
 * Management)](https://cloud.google.com/pubsub/access_control) allows you to
 * set permissions on invidual resources and offers a wider range of roles:
 * editor, owner, publisher, subscriber, and viewer. This gives you greater
 * flexibility and allows you to set more fine-grained access control.
 *
 * For example:
 *   * Grant access on a per-topic or per-subscription basis, rather than for
 *     the whole Cloud project.
 *   * Grant access with limited capabilities, such as to only publish messages
 *     to a topic, or to only to consume messages from a subscription, but not
 *     to delete the topic or subscription.
 *
 *
 * *The IAM access control features described in this document are Beta,
 * including the API methods to get and set IAM policies, and to test IAM
 * permissions. Cloud Pub/Sub's use of IAM features is not covered by any
 * SLA or deprecation policy, and may be subject to backward-incompatible
 * changes.*
 *
 * @class
 * @param {PubSub} pubsub PubSub Object.
 * @param {string} id The name of the topic or subscription.
 *
 * @see [Access Control Overview]{@link https://cloud.google.com/pubsub/access_control}
 * @see [What is Cloud IAM?]{@link https://cloud.google.com/iam/}
 *
 * @example
 * const {PubSub} = require('@google-cloud/pubsub');
 * const pubsub = new PubSub();
 *
 * const topic = pubsub.topic('my-topic');
 * // topic.iam
 *
 * const subscription = pubsub.subscription('my-subscription');
 * // subscription.iam
 */
export declare class IAM {
    Promise?: PromiseConstructor;
    pubsub: PubSub;
    request: typeof PubSub.prototype.request;
    id: string;
    constructor(pubsub: PubSub, id: string);
    getPolicy(gaxOpts?: CallOptions): Promise<GetPolicyResponse>;
    getPolicy(callback: GetPolicyCallback): void;
    getPolicy(gaxOpts: CallOptions, callback: GetPolicyCallback): void;
    setPolicy(policy: Policy, gaxOpts?: CallOptions): Promise<SetPolicyResponse>;
    setPolicy(policy: Policy, gaxOpts: CallOptions, callback: SetPolicyCallback): void;
    setPolicy(policy: Policy, callback: SetPolicyCallback): void;
    testPermissions(permissions: string | string[], gaxOpts?: CallOptions): Promise<TestIamPermissionsResponse>;
    testPermissions(permissions: string | string[], gaxOpts: CallOptions, callback: TestIamPermissionsCallback): void;
    testPermissions(permissions: string | string[], callback: TestIamPermissionsCallback): void;
}
