---
title: Two-factor authentication
subtitle: Enable two-factor authentication for the nf-core GitHub organisation
shortTitle: Two-factor authentication
---

The nf-core GitHub organisation requires two-factor authentication (2FA) for all members and outside collaborators.
2FA adds a second verification step when you sign in to GitHub, which protects nf-core repositories if your password is compromised.

## Who needs to act

What you need to do depends on your relationship with the nf-core GitHub organisation:

| Your situation                                             | What you need to do                                                                                                 |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Not a member of the nf-core GitHub organisation            | Nothing.                                                                                                            |
| Member of the nf-core GitHub organisation with 2FA enabled | Nothing.                                                                                                            |
| Member of the nf-core GitHub organisation without 2FA      | [Enable 2FA](#enable-two-factor-authentication).                                                                    |
| Outside collaborator on an nf-core repository without 2FA  | [Enable 2FA](#enable-two-factor-authentication) and [join the organisation](#join-the-nf-core-github-organisation). |

:::tip
To check whether you have 2FA enabled, go to **Settings** > **Password and authentication** on GitHub.
:::

## Enable two-factor authentication

GitHub supports several 2FA methods, including authenticator apps, security keys, passkeys, and GitHub Mobile.
Choose the method you prefer and follow the GitHub instructions to [configure two-factor authentication](https://docs.github.com/en/authentication/securing-your-account-with-two-factor-authentication-2fa/configuring-two-factor-authentication).

## Join the nf-core GitHub organisation

Outside collaborators are people with access to an nf-core repository who are not members of the nf-core GitHub organisation.
If you are an outside collaborator, request an invitation to the organisation in the [#github-invitations](https://nfcore.slack.com/channels/github-invitations) Slack channel.
See [Join nf-core](https://nf-co.re/join#github) for more information.

## What happens if you do not enable 2FA

The consequences of not enabling 2FA differ for members and outside collaborators:

- **Members** stay in the nf-core GitHub organisation, but cannot access organisation resources until they enable 2FA.
  This means you can no longer commit to, review, or merge pull requests in nf-core repositories.
- **Outside collaborators** are removed from the organisation and lose access to all nf-core repositories.

See [Requiring two-factor authentication in your organization](https://docs.github.com/en/organizations/keeping-your-organization-secure/managing-two-factor-authentication-for-your-organization/requiring-two-factor-authentication-in-your-organization) for details on how GitHub applies this requirement.

:::note
Anyone can still open a pull request from a fork without being a member of the nf-core GitHub organisation.
See [How to contribute](/docs/contributing/how-to-contribute) for more information.
:::
