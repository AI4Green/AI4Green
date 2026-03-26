import { useProfile } from "api";
import Cookies from "js-cookie";
import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";
import { useTranslation } from "react-i18next";

const UserContext = createContext({});

export const useUser = () => useContext(UserContext);

const getCookieProfile = () => {
  const yum = Cookies.get(".AI4Green4Students.Profile");
  return yum ? JSON.parse(yum) : null;
};

/**
 * Checks User Status on app load,
 * and provides methods to sign a user in and out
 * in response to app events (e.g. Login/Logout)
 */
export const UserProvider = ({ children }) => {
  const { i18n } = useTranslation();
  // const [user, setUser] = useState(getCookieProfile());

  const [user, setUser] = useState({
    fullName: "Guest User",
    email: "guest@example.com",
    uiCulture: "en",
  });

  const { data: profile, mutate } = useProfile();

  useEffect(() => {
    setUser(profile);
  }, [profile]);

  useEffect(() => {
    user && i18n.changeLanguage(user.uiCulture);
  }, [i18n, user]);

  const signOut = useCallback(() => setUser(null), []);
  const updateProfile = useCallback(() => mutate(), [mutate]);

  const context = useMemo(
    () => ({ user, signIn: setUser, signOut, updateProfile }),
    [user, setUser, signOut, updateProfile],
  );

  return (
    <UserContext.Provider value={context}>{children}</UserContext.Provider>
  );
};
