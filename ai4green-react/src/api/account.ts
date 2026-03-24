export const getAccountApi = ({ api }) => ({
  /**
   * Register a new User account
   * @param {*} body
   */
  register: (values) =>
    api.post("account/register", {
      json: values,
    }),

  /**
   * Login a user, given their username and password
   * @param {*} values
   */
  login: (values) =>
    api.post("account/login", {
      json: values,
    }),

  /**
   * Logout the current user
   */
  logout: () => api.post("account/logout"),

  /**
   * Try to confirm a User account
   * @param {*} userId The User ID
   * @param {*} token The previous Account Confirmation token
   * @returns
   */
  confirm: (userId, token) =>
    api.post("account/confirm", {
      json: { userId, token },
    }),

  /**
   * Resend an account confirmation email
   * @param {*} userIdOrEmail The User ID or Email Address
   * @returns
   */
  resendConfirm: (userIdOrEmail) =>
    api.post("account/confirm/resend", {
      json: userIdOrEmail,
    }),

  /**
   * Try to confirm a email change
   * @param {*} userId The User ID
   * @param {*} token The email change token
   * @param {*} newEmail The new email
   * @returns
   */
  confirmEmailChange: (userId, token, newEmail) =>
    api.post("account/email/confirm-change", {
      json: {
        credentials: { userId, token },
        data: { newEmail },
      },
    }),

  /**
   * Request a password reset email
   * @param {*} userIdOrEmail The User ID or Email Address
   * @returns
   */
  requestPasswordReset: (userIdOrEmail) =>
    api.post("account/password/request-reset", {
      json: userIdOrEmail,
    }),

  /**
   * Reset a User's password, using a valid token
   * @param {*} userId User ID to reset password for
   * @param {*} token System issued password reset token
   * @param {*} password the new password
   * @param {*} passwordConfirm confirm the new password
   * @returns
   */
  resetPassword: (userId, token, password, passwordConfirm) =>
    api.post("account/password/reset", {
      json: {
        credentials: { userId, token },
        data: { password, passwordConfirm },
      },
    }),

  /**
   * Activate Users's account, using a valid token
   * @param {*} userId User ID to activate account for
   * @param {*} token System issued account activation token
   * @param {*} password the password
   * @param {*} fullName User full name
   * @returns
   */
  activateAccount: ({ userId, token, password, fullName }) =>
    api.put("account/activate", {
      json: {
        credentials: { userId, token },
        data: { password, fullName },
      },
    }),

  /**
   * Login with OIDC provider.
   * @param {*} redirectUri The redirect URI.
   * @param {*} idp The OIDC provider. e.g. "Microsoft", must match the alias set in Keycloak.
   * @returns
   */
  oidcLogin: ({ redirectUri, idp }) => {
    const url = new URL("/api/account/oidc-login", window.location.origin);
    url.searchParams.set("redirectUri", redirectUri);
    url.searchParams.set("idp", idp);

    window.location.href = url.toString();
  },
});
