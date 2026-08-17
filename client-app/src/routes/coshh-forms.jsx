import { NotFound } from "pages/error";
import { ProjectTypeCanvas } from "pages/project-type";
import { Route, Routes } from "react-router-dom";
import { CoshhFormsList } from "pages/project-type";

export const CoshhForms = () => {
  return (
    <Routes>
      <Route path="/" element={<CoshhFormsList />} />
      <Route path="/:projectTypeId">
        <Route index element={<ProjectTypeCanvas />} />
        <Route path="sections/:sectionId" element={<ProjectTypeCanvas />} />
        <Route
          path="section-types/:sectionTypeId/sections/:sectionId"
          element={<ProjectTypeCanvas />}
        />
      </Route>
      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};
