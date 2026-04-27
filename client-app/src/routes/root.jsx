import { DefaultLayout } from "layouts/default";
import { ContentPage } from "pages/content";
import { NotFound } from "pages/error";
import { Home } from "pages/home";
import { Route, Routes } from "react-router-dom";
import { ProjectType } from "./project-type";
import { Project } from "./project";

export const Root = () => {
  return (
    <Routes>
      {/* 1. Home page usually stands alone */}
      <Route path="/" element={<Home />} />

      {/* 2. Routes that share the DefaultLayout */}
      <Route path="/" element={<DefaultLayout />}>
        <Route
          path="greenchemistry"
          element={<ContentPage contentKey="greenchemistry" />}
        />
        <Route path="about" element={<ContentPage contentKey="about" />} />
        <Route
          path="documentation"
          element={<ContentPage contentKey="documentation" />}
        />

        <Route path="project-types/*" element={<ProjectType />} />
        <Route path="projects/*" element={<Project />} />
      </Route>

      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};
