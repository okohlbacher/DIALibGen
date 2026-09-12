# Every file that carries this project's version must agree with
# project(DIALibGen VERSION ...) in CMakeLists.txt:
#
#   gui/package.json               the frontend package
#   gui/package-lock.json          npm's two copies of package.json's version
#   gui/src-tauri/tauri.conf.json  stamps the desktop bundle (.app / .dmg / .msi)
#   gui/src-tauri/Cargo.toml       the Tauri backend crate
#   gui/src-tauri/Cargo.lock       cargo's copy of Cargo.toml's version
#
# Nothing at run time reads most of these, so they drift in silence: the
# reference implementation this was taken from shipped three releases with the
# Cargo.toml version stuck two minors behind, and one release announcing
# OpenMS's version instead of its own. Checking at CONFIGURE time means the
# first cmake after a bump fails -- on a developer's machine and on every CI
# leg, all of which configure.
#
# Versions are plain MAJOR.MINOR.PATCH everywhere: project(VERSION) accepts
# nothing else, and the comparison is textual, so "0.2.0.0" or "0.2.0-rc1" do
# not pass as 0.2.0.
#
# Skipped entirely when gui/ is absent (a CLI-only source package). A single
# missing file inside a present gui/ is an error, not a skip.
#
# To bump: edit CMakeLists.txt, gui/package.json, gui/src-tauri/tauri.conf.json
# and gui/src-tauri/Cargo.toml, then refresh the locks with
# `npm install --package-lock-only` in gui/ and
# `cargo update -p dialibgen-gui --offline` in gui/src-tauri/.

function(_dlg_check file want values)
  list(LENGTH values n)
  if(n EQUAL 0)
    message(FATAL_ERROR "version guard: no version found in ${file}")
  endif()
  foreach(v IN LISTS values)
    if(NOT v STREQUAL want)
      message(FATAL_ERROR
        "version guard: ${file} says ${v} but CMakeLists.txt says ${want}. "
        "See cmake/version-guard.cmake for what a bump has to touch.")
    endif()
  endforeach()
endfunction()

function(dialibgen_version_guard want root)
  if(NOT IS_DIRECTORY "${root}/gui")
    message(STATUS "DIALibGen: no gui/ in this source tree -- version guard skipped")
    return()
  endif()

  set(files
    gui/package.json
    gui/package-lock.json
    gui/src-tauri/tauri.conf.json
    gui/src-tauri/Cargo.toml
    gui/src-tauri/Cargo.lock)
  foreach(f IN LISTS files)
    if(NOT EXISTS "${root}/${f}")
      message(FATAL_ERROR "version guard: ${f} is missing")
    endif()
    # Re-run when one of THEM changes, not only when CMakeLists.txt does.
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${root}/${f}")
  endforeach()

  # --- the JSON ones: a real parser, so formatting cannot fool it -----------
  foreach(f IN ITEMS gui/package.json gui/package-lock.json gui/src-tauri/tauri.conf.json)
    file(READ "${root}/${f}" text)
    set(found "")
    string(JSON v ERROR_VARIABLE err GET "${text}" version)
    if(NOT err)
      list(APPEND found "${v}")
    endif()
    if(f STREQUAL "gui/package-lock.json")
      # npm writes the root version a second time, under packages[""].
      string(JSON v ERROR_VARIABLE err GET "${text}" packages "" version)
      if(NOT err)
        list(APPEND found "${v}")
      endif()
      list(LENGTH found n)
      if(NOT n EQUAL 2)
        message(FATAL_ERROR "version guard: ${f} carries ${n} versions, expected 2")
      endif()
    endif()
    _dlg_check("${f}" "${want}" "${found}")
  endforeach()

  # --- the cargo ones ------------------------------------------------------
  # Cargo writes `name` immediately above `version` in both files, so anchoring
  # on our own crate name picks our block out of Cargo.lock's hundreds without
  # a TOML parser. A CRLF checkout is normalised first.
  foreach(f IN ITEMS gui/src-tauri/Cargo.toml gui/src-tauri/Cargo.lock)
    file(READ "${root}/${f}" text)
    string(REPLACE "\r" "" text "${text}")
    if(NOT text MATCHES "name[ \t]*=[ \t]*\"dialibgen-gui\"[ \t]*\n+version[ \t]*=[ \t]*\"([^\"]*)\"")
      message(FATAL_ERROR
        "version guard: no `name = \"dialibgen-gui\"` followed by a version in ${f}")
    endif()
    _dlg_check("${f}" "${want}" "${CMAKE_MATCH_1}")
  endforeach()

  list(JOIN files ", " files)
  message(STATUS "DIALibGen: version ${want} agrees across CMakeLists.txt, ${files}")
endfunction()

dialibgen_version_guard("${PROJECT_VERSION}" "${CMAKE_CURRENT_SOURCE_DIR}")
