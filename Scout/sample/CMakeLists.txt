cmake_minimum_required(VERSION 3.5)
project(ScoutMinimalSample)

add_executable(${PROJECT_NAME} 
    main.cpp
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    Scout
)
