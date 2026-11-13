import { DefaultLayout } from "layouts/default";
import { ContentPage } from "pages/content";
import { NotFound } from "pages/error";
import { Home } from "pages/home";
import { Route, Routes } from "react-router-dom";
import { CoshhForms } from "./coshh-forms";
import { COSHH } from "./coshh.jsx";
import InventoryPage from "components/inventory/inventory-page.jsx";
import LocationsPage from "components/inventory/locations-page.jsx";

export const Root = () => {
  return (
    <Routes>
      <Route path="/" element={<Home />} />

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

        <Route path="coshh-form-management/*" element={<CoshhForms />} />
        <Route path="my-coshh-forms/*" element={<COSHH />} />

        <Route path="inventory/:workgroupName" element={<InventoryPage />} />
        <Route
          path="inventory/:workgroupName/locations"
          element={<LocationsPage />}
        />
      </Route>

      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};
