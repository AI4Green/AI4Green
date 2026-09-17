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

  const initialUser = getCookieProfile();

  const [user, setUser] = useState(initialUser);
  const [isLoading, setIsLoading] = useState(!initialUser);

  const { data: profile, mutate } = useProfile();

  useEffect(() => {
    if (profile === undefined) return;

    // API says logged in
    if (profile) {
      setUser(profile);

      Cookies.set(".AI4Green4Students.Profile", JSON.stringify(profile), {
        expires: 1,
      });
    }
    // API says NOT logged in
    else {
      setUser(null);
      Cookies.remove(".AI4Green4Students.Profile");
    }

    setIsLoading(false);
  }, [profile]);

  useEffect(() => {
    if (user?.uiCulture) {
      i18n.changeLanguage(user.uiCulture);
    }
  }, [i18n, user]);

  const signOut = useCallback(() => {
    setUser(null);
    Cookies.remove(".AI4Green4Students.Profile");
  }, []);

  const updateProfile = useCallback(() => mutate(), [mutate]);

  const value = useMemo(
    () => ({
      user,
      isLoading,
      signIn: setUser,
      signOut,
      updateProfile,
    }),
    [user, isLoading, signOut, updateProfile],
  );

  return <UserContext.Provider value={value}>{children}</UserContext.Provider>;
};
